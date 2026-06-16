import ast
import concurrent
import json
import logging
import os
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from fdr import prob_fdr
import pandas as pd
import yaml


sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import utils, global_data_process as gdp, constants as const


class COLOC2QTLLOCI:

    TOOL_NAME = 'coloc'

    def __init__(self):
        logging.info('init MatchingProcessor for loci')

    def run(self, processor):

        self.processor = processor
        self.config_holder = self.processor.config_holder
        self.global_config = self.processor.global_config
        self.root_work_dir = self.processor.root_work_dir
        self.output_preprocessed_dir = self.processor.output_preprocessed_dir
        self.output_processed_dir = self.processor.output_processed_dir


        self.tool_parent_dir = self.processor.tool_parent_dir
        self.rank_dir = self.processor.rank_dir
        # qtl1
        self.qtl1_type = self.processor.qtl1_type
        self.qtl1_biological_context = self.processor.qtl1_biological_context
        self.qtl1_col_dict = self.processor.qtl1_col_dict
        self.qtl1_preprocesed_dir = self.processor.qtl1_preprocesed_dir
        self.qtl1_grouped_dir = self.processor.qtl1_grouped_dir
        self.qtl1_output_report = self.processor.qtl1_output_report
        self.qtl1_LD_window = self.processor.qtl1_LD_window

        qtl1_type_dict = {self.qtl1_col_dict['chrom']: 'category',
                          self.qtl1_col_dict['position']: 'Int64',
                          self.qtl1_col_dict['alt']: 'category',
                          self.qtl1_col_dict['ref']: 'category',
                          self.qtl1_col_dict['phenotype_id']: 'category'
                          }
        self.qtl1_sample_size = self.processor.global_config['input']['qtl1'].get('sample_size', 948)
        self.qtl1_threshold = self.processor.config_holder.qtl1_p_threshold

        # qtl2
        self.qtl2_type = self.processor.qtl2_type
        self.qtl2_biological_context = self.processor.qtl2_biological_context
        self.qtl2_preprocesed_dir = self.processor.qtl2_preprocesed_dir
        self.qtl2_grouped_dir = self.processor.qtl2_grouped_dir
        self.qtl2_output_report = self.processor.qtl2_output_report
        self.qtl2_LD_window = self.processor.qtl2_LD_window
        self.qtl2_col_dict = self.processor.qtl2_col_dict

        qtl2_type_dict = {self.qtl2_col_dict['chrom']: 'category',
                          self.qtl2_col_dict['position']: 'Int64',
                          self.qtl2_col_dict['alt']: 'category',
                          self.qtl2_col_dict['ref']: 'category',
                          self.qtl2_col_dict['phenotype_id']: 'category'
                          }
        self.qtl2_sample_size = self.processor.global_config['input']['qtl2'].get('sample_size', 948)
        self.qtl2_threshold = self.processor.config_holder.qtl2_p_threshold

        # vcf input output path
        self.ref_vcf_dir = self.processor.ref_vcf_dir
        self.rsidvcf = self.processor.rsidvcf
        self.min_matching_number = self.processor.min_matching_number
        print(f"globalprocessorself.rsidvcf: {self.rsidvcf}")

        self.population = self.config_holder.population
        self.shell_command_plink_execute = 'plink --silent --vcf {} --r2 --matrix --mac 1 --write-snplist --out {}'

        # tools parameter config path
        self.tools_config_file = self.processor.tools_config_file
        self.var_id_col_name = self.processor.VAR_ID_COL_NAME

        self.coloc_base_dir = self.processor.coloc_base_dir

        qtl1_summary_df = pd.read_csv(self.qtl1_output_report, sep=const.column_spliter,
                                      dtype={self.qtl1_col_dict['chrom']: 'category'})
        qtl2_summary_df = pd.read_csv(self.qtl2_output_report, sep=const.column_spliter,
                                      dtype={self.qtl2_col_dict['chrom']: 'category'})

        if not Path(os.path.join(self.tool_parent_dir, f'{self.TOOL_NAME}_processed_done.txt')).exists():
            logging.info(f"COLOC qtl_qtl_loci_processing")
            start_time = datetime.now()
            Path(self.output_processed_dir).mkdir(parents=True, exist_ok=True)
            logging.info(f'run_coloc start at: {start_time}')


            ## coloc/input
            self.coloc_dir_input = os.path.join(self.coloc_base_dir, 'input')   # combined_ld/coloc/input
            Path(self.coloc_dir_input).mkdir(parents=True, exist_ok=True)

            self.start_process(qtl1_summary_df, self.qtl1_col_dict, qtl2_summary_df, self.qtl2_col_dict)

            with open(os.path.join(self.tool_parent_dir, f'{self.TOOL_NAME}_processed_done.txt'), 'w') as f:
                f.write("DONE\n")


    def start_process(self, qtl1_summary_df, qtl1_col_dict, qtl2_summary_df, qtl2_col_dict):
        logging.info(f"start_process")
        total_len = len(qtl2_summary_df)
        for qtl1_ix, qtl1_row in qtl1_summary_df.iterrows():
            qtl1_chrom = str(qtl1_row.loc['chrom'])
            for qtl2_ix, qtl2_row in qtl2_summary_df.iterrows():
                qtl2_chrom = str(qtl2_row.loc['chrom'])
                if qtl1_chrom != qtl2_chrom:
                    logging.info(f"qtl chrom not in match")
                    continue
                chrom = qtl1_chrom
                print(f"processing {chrom}")
                
                # logging.info(f"qtl_gene_file: {qtl_gene_file}")
                qtl1_phenotype_id = qtl1_row.loc['pheno_file'].strip('.tsv.gz')
                qtl2_phenotype_id = qtl2_row.loc['pheno_file'].strip('.tsv.gz')
                # logging.info(f"phenotype_id: {phenotype_id}")

                qtl1_significant_positions = ast.literal_eval(qtl1_row.loc['positions'])
                qtl2_significant_positions = ast.literal_eval(qtl2_row.loc['positions'])

                if len(set(qtl1_significant_positions) & set(qtl2_significant_positions)) < self.min_matching_number:
                    logging.info(f"coloc not enough SNPs in {self.qtl1_type} {qtl1_phenotype_id} and {self.qtl2_type} {qtl2_phenotype_id} ")
                    continue
                print(f"checkpoint 4: {qtl1_phenotype_id}_{qtl2_phenotype_id}")
                qtl1_pheno_file = os.path.join(self.qtl1_grouped_dir, chrom, f"{qtl1_phenotype_id}.tsv.gz")
                qtl2_pheno_file = os.path.join(self.qtl2_grouped_dir, chrom, f"{qtl2_phenotype_id}.tsv.gz")
                self.process_pheno(qtl1_pheno_file, qtl1_col_dict, qtl1_phenotype_id,
                                   qtl2_pheno_file, qtl2_col_dict, qtl2_phenotype_id,
                                    self.var_id_col_name, chrom, 
                                    self.qtl1_sample_size, self.qtl2_sample_size,
                                    self.min_matching_number, 
                                    self.qtl1_threshold, self.qtl2_threshold)


    def set_gwas_range_files(self, gwas_cluster_output_dir):
        gwas_range_files = {}
        for gwas_range_file in os.listdir(gwas_cluster_output_dir):
            part_list = utils.get_file_name(gwas_range_file).split('_')
            # if len(part_list) < 2 or 'chr' not in part_list[1]:
            #     continue
            if len(part_list) < 2 or 'chr' not in part_list[0]:
                continue
            chrom = part_list[0].strip('chr')
            range_files = gwas_range_files.get(chrom, [])
            range_files.append(os.path.join(gwas_cluster_output_dir, gwas_range_file))
            gwas_range_files[chrom] = range_files
        
        return gwas_range_files



    def process_pheno(self, qtl1_pheno_file, qtl1_col_dict, qtl1_phenotype_id,
                      qtl2_pheno_file, qtl2_col_dict, qtl2_phenotype_id,
                      var_id_col_name, chrom, 
                      qtl1_sample_size, qtl2_sample_size,
                      min_matching_number, 
                      qtl1_threshold, qtl2_threshold):
        logging.info(f'COLOC process_pheno {self.qtl1_type} {qtl1_phenotype_id} and {self.qtl2_type} {qtl2_phenotype_id}')
        # print(f'COLOC process_pheno')

        qtl1_trait_df = pd.read_table(qtl1_pheno_file, sep=const.column_spliter,
                                          usecols=[
                                              var_id_col_name,
                                              qtl1_col_dict['chrom'],
                                              qtl1_col_dict['position'],
                                              qtl1_col_dict['alt'],
                                              qtl1_col_dict['ref'],
                                              qtl1_col_dict['beta'],
                                              qtl1_col_dict['se'],
                                              qtl1_col_dict['pvalue'],
                                              qtl1_col_dict['snp'],
                                              qtl1_col_dict['phenotype_id'],
                                              qtl1_col_dict['maf'], 
                                              ],
                                          dtype=qtl1_col_dict)
        if len(qtl1_trait_df) <= 1:
            logging.info(f'no qtl1_trait_df')
            return

        if not Path(qtl2_pheno_file).exists():
            logging.info(f"ProcessingPheno:: {qtl2_pheno_file} does not exists!!")
            return
        print(f"checkpoint 3: {qtl2_phenotype_id}")
        qtl2_trait_df = pd.read_table(qtl2_pheno_file, sep=const.column_spliter,
                                      usecols=[
                                          var_id_col_name,
                                          qtl2_col_dict['chrom'],
                                          qtl2_col_dict['position'],
                                          qtl2_col_dict['alt'],
                                          qtl2_col_dict['ref'],
                                          qtl2_col_dict['beta'],
                                          qtl2_col_dict['se'],
                                          qtl2_col_dict['pvalue'],
                                          qtl2_col_dict['phenotype_id'],
                                          qtl2_col_dict['maf'],
                                          qtl2_col_dict['snp']],
                                      dtype=qtl2_col_dict)

        qtl2_trait_df.drop(
            index=qtl2_trait_df[~qtl2_trait_df[var_id_col_name].isin(qtl1_trait_df[var_id_col_name])].index,
            inplace=True)
        # print(f"qtl2_trait_df[var_id_col_name]: {qtl2_trait_df[var_id_col_name]}")
        # print(f"candidate_gwas_df[var_id_col_name]: {candidate_gwas_df[var_id_col_name]}")
        
        if len(qtl2_trait_df) < min_matching_number:
            logging.info(f'Skip, qtl no more than {min_matching_number}') # 有问题
            return
        # print(f"len(qtl_trait_df): {len(qtl_trait_df)}")
        utils.drop_non_intersect_rows(qtl2_trait_df, var_id_col_name, qtl1_trait_df, var_id_col_name)
        if len(qtl1_trait_df) < min_matching_number:
            logging.info(f'Skip, gwas no more than {min_matching_number}')
            return
        ## coloc/input

        coloc_input_qtl1_path = os.path.join(self.coloc_dir_input, f'qtl1_{self.qtl1_type}_{chrom}_{qtl1_phenotype_id}.tsv.gz')
        coloc_input_qtl2_path = os.path.join(self.coloc_dir_input, f'qtl2_{self.qtl2_type}_{chrom}_{qtl2_phenotype_id}.tsv.gz')
        print(f'coloc_qtl1_input_path: {coloc_input_qtl1_path}')
        print(f'coloc_qtl2_input_path: {coloc_input_qtl2_path}')



        # Adjust qtl beta sign according to gwas allele order
        qtl1_trait_df.reset_index(drop=True, inplace=True)
        qtl2_trait_df.reset_index(drop=True, inplace=True)
        # utils.adjust_allele_order(candidate_gwas_df,
        #                           gwas_col_dict['effect_allele'],
        #                           gwas_col_dict['other_allele'],
        #                           gwas_col_dict['chrom'],
        #                           gwas_col_dict['position'],
        #                           qtl_trait_df,
        #                           ref_df_chrom_col_name=qtl_col_dict['chrom'],
        #                           ref_df_pos_col_name=qtl_col_dict['position'],
        #                           ref_df_alt_allele_col_name=qtl_col_dict['alt'],
        #                           ref_df_ref_allele_col_name=qtl_col_dict['ref'],
        #                           gbeta_col_name=gwas_col_dict['beta'])


        # Now gwas/qtl/vcf have the same num of rows, write candidate data to file
        # Reverse GWAS&qtl column mapping key-value and pass to R so that dataframe in R has fixed column names
        qtl2_trait_df_coloc = qtl2_trait_df[[var_id_col_name,
                                             qtl2_col_dict['chrom'],
                                             qtl2_col_dict['position'], 
                                             qtl2_col_dict['alt'],
                                             qtl2_col_dict['ref'],
                                             qtl2_col_dict['beta'],
                                             qtl2_col_dict['se'],
                                             qtl2_col_dict['pvalue'],
                                             qtl2_col_dict['phenotype_id'],
                                             qtl2_col_dict['maf']]]
        qtl1_trait_df_coloc = qtl1_trait_df[[var_id_col_name,
                                             qtl1_col_dict['chrom'],
                                             qtl1_col_dict['position'], 
                                             qtl1_col_dict['alt'],
                                             qtl1_col_dict['ref'],
                                             qtl1_col_dict['beta'],
                                             qtl1_col_dict['se'],
                                             qtl1_col_dict['pvalue'],
                                             qtl1_col_dict['phenotype_id'],
                                             qtl1_col_dict['maf']]]

        if ('varbeta' not in qtl1_col_dict.keys() or qtl1_col_dict.get('varbeta') is None) and (
                'se' in qtl1_col_dict.keys() and qtl1_col_dict.get('se') is not None):
            qtl1_trait_df_coloc['varbeta'] = qtl1_trait_df_coloc[qtl1_col_dict['se']] ** 2
        qtl1_trait_df_coloc.rename({v: k for k, v in qtl1_col_dict.items()}, axis='columns', inplace=True)
        if ('varbeta' not in qtl2_trait_df_coloc.keys() or qtl2_col_dict.get('varbeta') is None) and (
                'se' in qtl2_col_dict.keys() and qtl2_col_dict.get('se') is not None):
            qtl2_trait_df_coloc['varbeta'] = qtl2_trait_df_coloc[qtl2_col_dict['se']] ** 2
        qtl2_trait_df_coloc.rename({v: k for k, v in qtl2_col_dict.items()}, axis='columns', inplace=True)


        if len(qtl1_trait_df_coloc[qtl1_trait_df_coloc['pvalue'] < qtl1_threshold]) <= 0:
            logging.info(f'{qtl1_phenotype_id}c oloc no sig qtl_trait {qtl1_threshold}')
            return
        if len(qtl2_trait_df_coloc[qtl2_trait_df_coloc['pvalue'] < qtl2_threshold]) <= 0:
            logging.info(f'{qtl2_phenotype_id} coloc no sig qtl_trait {qtl2_threshold}')
            return

        print(f"successfully output coloc input files for {self.qtl1_type} {qtl1_phenotype_id} and {self.qtl2_type} {qtl2_phenotype_id}")
        qtl1_trait_df_coloc.to_csv(\
                            coloc_input_qtl1_path, sep=const.output_spliter, header=True, index=False)
        qtl2_trait_df_coloc.to_csv(\
                            coloc_input_qtl2_path, sep=const.output_spliter, header=True, index=False)

    
    def __convert_positions_str_to_list(self, positions_str):

        if isinstance(positions_str, str):
            return json.loads(positions_str)
        elif isinstance(positions_str, list):
            return positions_str
        else:
            return []

    def __get_cluster_significant_snps_dict(self, gwas_summary_df):

        cluster_snps_dict = {}
        for _, row in gwas_summary_df.iterrows():
            key = row.loc['lead_SNP'] + '_chr' + str(row.loc['chrom']) + '_' + str(row.loc['range_start']) + '_' + str(row.loc['range_end'])
            cluster_snps_dict[key] = self.__convert_positions_str_to_list(row.loc['positions'])
        return cluster_snps_dict
