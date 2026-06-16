import concurrent
import datetime
import logging
import os
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import ast

import pandas as pd

sys.path.append(os.path.abspath(os.path.dirname(Path(__file__).resolve())))
import utils
import constants as const

'''
预处理的GWAS文件路径:
{working_dir}/preprocessed/gwas/{trait}
预处理的eQTL文件路径:
{working_dir}/preprocessed/eqtl/{tissue}

工具working dir:
{working_dir}/processed/{trait}_{tissue}/{tool}
{working_dir}/processed/{trait}_{tissue}/{tool}
'''


class Processor:
    # column to identify a unique record in gwas & eqtl files, in chromosome_position format
    VAR_ID_COL_NAME = 'var_id_'
    def __init__(self, cfg_holder=None):
        logging.info('init Processor')
        self.config_holder = cfg_holder
        self.global_config = self.config_holder.global_config
        self.root_work_dir = self.config_holder.root_work_dir
        self.tools = self.global_config['tools']
        print(f"global_config: {self.global_config}")
        self.output_preprocessed_dir = self.config_holder.output_preprocessed_dir
        self.output_processed_dir = self.config_holder.output_processed_dir
        self.min_matching_number = self.config_holder.min_matching_number
        Path(self.output_preprocessed_dir).mkdir(exist_ok=True, parents=True)
        Path(self.output_processed_dir).mkdir(exist_ok=True, parents=True)     
        self.study_dir = self.config_holder.study_dir

        self.genomic_window = self.config_holder.genomic_window
        self.window_size = self.config_holder.window_size
        self.LD_r2_filter = self.config_holder.LD_r2_filter
        self.LD_additional_expansion = self.config_holder.LD_additional_expansion
        self.population = self.config_holder.population
        self.genecode_file = self.config_holder.genecode_file
        self.report_path = self.config_holder.report_path

        self.ref_vcf_dir = self.config_holder.ref_vcf_dir
        self.rsidvcf = self.config_holder.rsidvcf
        self.ld_block_vcf_dir = self.config_holder.ld_block_vcf_dir
        self.parallel = self.config_holder.parallel
        self.results = self.config_holder.results

        self.keep_intermediate = self.config_holder.keep_intermediate
        self.tools_config_file = self.config_holder.tools_config_file
        self.whether_fixed_window = self.config_holder.whether_fixed_window

        self.cCREs_file = self.config_holder.cCREs_file
        ########################################################################
        ##                            GWAS setting                            ##
        ########################################################################



        ###### qtl1
        self.tool_parent_dir = self.config_holder.tool_parent_dir
        Path(self.tool_parent_dir).mkdir(exist_ok=True, parents=True)


        self.qtl1_type = self.config_holder.qtl1_type
        self.qtl1_file = self.config_holder.qtl1_file
        self._qtl1_type = self.config_holder._qtl1_type
        self.qtl1_biological_context = self.config_holder.qtl1_biological_context

        self.qtl1_grouped_dir = self.config_holder.qtl1_grouped_dir
        self.qtl1_output_report = self.config_holder.qtl1_output_report
        self.qtl1_LD_window = self.config_holder.qtl1_LD_window
        self.qtl1_preprocesed_dir = self.config_holder.qtl1_preprocesed_dir
        self.qtl1_col_dict = self.config_holder.qtl1_col_dict
        self.qtl1_p_threshold = self.config_holder.qtl1_p_threshold
        self.qtl1_sep = self.config_holder.qtl1_sep
        self.qtl1_sample_size = self.config_holder.qtl1_sample_size

        ###### qtl2

        self.qtl2_type = self.config_holder.qtl2_type
        self.qtl2_file = self.config_holder.qtl2_file
        self._qtl2_type = self.config_holder._qtl2_type
        self.qtl2_biological_context = self.config_holder.qtl2_biological_context

        self.qtl2_grouped_dir = self.config_holder.qtl2_grouped_dir
        self.qtl2_output_report = self.config_holder.qtl2_output_report
        self.qtl2_LD_window = self.config_holder.qtl2_LD_window
        self.qtl2_preprocesed_dir = self.config_holder.qtl2_preprocesed_dir
        self.qtl2_col_dict = self.config_holder.qtl2_col_dict
        self.qtl2_p_threshold = self.config_holder.qtl2_p_threshold
        self.qtl2_sep = self.config_holder.qtl2_sep
        self.qtl2_sample_size = self.config_holder.qtl2_sample_size


        self.coloc_base_dir = self.config_holder.coloc_base_dir
        self.susie_base_dir = self.config_holder.susie_base_dir
        self.rank_dir = self.config_holder.rank_dir

        self.shell_command_plink_execute = 'plink --silent --vcf {} --r2 gz --ld-snp {} --ld-window 99999 --ld-window-kb {} --ld-window-r2 {} --out {}'
        self.shell_command_plink_ldblock_execute = 'plink --silent --vcf {} --r2 gz --ld-window 999999 --ld-window-r2 {} --out {}'


    def __merge_alt_ref(self, gwas_df, chrom, population):
        logging.info(f'Merging alt ref for chrom {chrom}')
        input_vcf = os.path.join(self.ref_vcf_dir, population, f'chr{chrom}.vcf.gz')
        vcf_df = pd.read_table(input_vcf, header=None, comment='#', usecols=[0, 1, 3, 4],
                               dtype={0: 'category', 
                                      1: 'Int64', 
                                    #   3: pd.CategoricalDtype(const.SNP_ALLELE),
                                    #   4: pd.CategoricalDtype(const.SNP_ALLELE)
                                      3: 'category',
                                      4: 'category'
                                      })
        vcf_df.columns = [self.gwas_col_dict['chrom'], self.gwas_col_dict['position'], 'ref_', 'alt_']
        vcf_df.dropna(subset=['ref_', 'alt_'], inplace=True)
        vcf_df.reset_index(drop=True, inplace=True)
        merged_vcf_col_suffix = f'_{chrom}_vcf'
        print(f"merged_vcf_col_suffix: {merged_vcf_col_suffix}")
        # print(f"__merge_alt_ref gwas_df: {gwas_df}")
        vcf_df['ref_'] = vcf_df['ref_'].str.upper()
        vcf_df['alt_'] = vcf_df['alt_'].str.upper()
        vcf_df.index = vcf_df[self.gwas_col_dict['chrom']].astype(str) + '_' + vcf_df[self.gwas_col_dict['position']].astype(str) + '_' + vcf_df['ref_'] + '_' + vcf_df['alt_']
        gwas_df.index = gwas_df[self.gwas_col_dict['variant_id']]
        gwas_df = pd.merge(left=gwas_df,
                           right=vcf_df,
                           left_index=True,
                           right_index=True,
                           how='left',
                           suffixes=(None, merged_vcf_col_suffix))
        del vcf_df
        gwas_df['ref_'].mask((gwas_df[self.gwas_col_dict['chrom']] == chrom) & gwas_df['ref_'].isna(),
                             gwas_df[f'ref_{merged_vcf_col_suffix}'], inplace=True)
        gwas_df['alt_'].mask((gwas_df[self.gwas_col_dict['chrom']] == chrom) & gwas_df['alt_'].isna(),
                             gwas_df[f'alt_{merged_vcf_col_suffix}'], inplace=True)
        gwas_df.drop(columns=[col for col in gwas_df.columns if col.endswith(merged_vcf_col_suffix)], inplace=True)
        logging.info(f'Merge alt ref for chrom {chrom} completed')

################################################################################
#                             1.3 Preprocess QTL                              #
################################################################################

    def __qtl_ld_interval(self, pval_filter_df, total_df, window_size, qtl_preprocesed_dir,
                          col_dict, chrom_list, lead_SNP_list, positions_list, 
                          cluster_pos_dict, count, p_threshold, outputpath, phenotype_id,
                          additional_expansion = 50000):

        population = self.population
        print(f"population: {population}")


        for _, row in pval_filter_df.iterrows(): 
            chrom = row.loc[col_dict['chrom']]
            pos = row.loc[col_dict['position']]
            chrom_cluster_pos_list = cluster_pos_dict.get(chrom, [])
            
            # print(f"chrom {chrom}, pos {pos}, chrom_cluster_pos_list {chrom_cluster_pos_list}")
            if pos in set(chrom_cluster_pos_list):
                print(f"position {pos} in chrom {chrom} already in cluster_pos_dict, continue")
                continue
            
            cluster_start = max(0, int(pos) - int(window_size))
            cluster_end = int(pos) + int(window_size)

            range_df = total_df[(total_df[col_dict['chrom']] == chrom) & (
                    cluster_start <= total_df[col_dict['position']]) & (
                                    total_df[col_dict['position']] <= cluster_end)]
            ############################# To speed up ###########################
            ## If the SNPs in range_df overlapped with cluster_pos_dict[chrom], regard this range_df
            if len(set(range_df[col_dict['position']]) & set(chrom_cluster_pos_list)) > 0:
                continue
            #####################################################################
            lead_SNP_id = row.loc[Processor.VAR_ID_COL_NAME]

            print(f"__qtl_ld_interval: {lead_SNP_id}") # chr1_13550_G_A

            vcf_output_dir = os.path.join(qtl_preprocesed_dir, 'vcf')
            # utils.tmp_delete_dir(vcf_output_dir)
            Path(vcf_output_dir).mkdir(parents=True, exist_ok=True)
            # output_vcf_name = f"{lead_SNP_id}.vcf"

            input_vcf = os.path.join(self.ref_vcf_dir, population, f'chr{chrom}.vcf.gz')

            ldblock_vcf = None
            ld_block_path = None
            if self.ld_block_vcf_dir:
                ldblock_vcf = utils.get_ref_ldblock_vcf_path(population, chrom, pos, self.ld_block_vcf_dir)
                ld_block_path = utils.get_ref_ldblock_ld_path(population, chrom, pos, self.ld_block_vcf_dir)
            if ldblock_vcf is not None:
                input_vcf = ldblock_vcf
            if lead_SNP_id == None:
                continue
            print(f"ld_block_path: {ld_block_path}")

            if ld_block_path and Path(ld_block_path).is_file():
                ld_path = ld_block_path
            else:
                output_ld_prefix = os.path.join(self.vcf_output_dir, lead_SNP_id)
                os.system(self.shell_command_plink_execute.format(input_vcf, 
                                                        lead_SNP_id, window_size/1000, self.LD_r2_filter,
                                                        output_ld_prefix))
                # os.system(self.shell_command_plink_ldblock_execute.format(input_vcf, self.LD_r2_filter, output_ld_block_prefix))
                ld_path = output_ld_prefix + '.ld.gz'

            if ld_path and Path(ld_path).is_file():
                print(f"DEBUG: ld_path: {ld_path}")
                ld_filter_df = pd.read_csv(ld_path, sep='\s+', usecols=['SNP_A','SNP_B','R2'])
                
                if (len(ld_filter_df) < 2) or (lead_SNP_id not in list(ld_filter_df['SNP_A'])) or (lead_SNP_id not in list(ld_filter_df['SNP_B'])):
                    logging.info(f"Calculating LD failure {lead_SNP_id}.")
                    range_df = range_df[(
                            max(0, pos - additional_expansion) <= range_df[col_dict['position']]) & (
                                            range_df[col_dict['position']] <= pos + additional_expansion)]

                else:
                    related_snps = pd.concat([
                        ld_filter_df.loc[ld_filter_df['SNP_A'] == lead_SNP_id, 'SNP_B'],
                        ld_filter_df.loc[ld_filter_df['SNP_B'] == lead_SNP_id, 'SNP_A'],
                    ]).to_frame('SNP')
                    related_snps = related_snps.drop_duplicates()
                    related_snps['pos'] = related_snps['SNP'].str.split('_').str[1].astype(int)

                    min_snp = related_snps.loc[related_snps['pos'].idxmin()]['SNP']
                    print(f"DEBUG: min_snp: {min_snp}")
                    max_snp = related_snps.loc[related_snps['pos'].idxmax()]['SNP']
                    print(f"DEBUG: max_snp: {max_snp}")

                    range_start = int(min_snp.split('_')[1])
                    range_end = int(max_snp.split('_')[1])
                    print(f"DEBUG: range: {range_start} - {range_end}")

                    range_df[col_dict['position']] = range_df[col_dict['position']].astype('int')
                    logging.info(f"DEBUG: print(range_df[col_dict['position']].dtype): {range_df[col_dict['position']].dtype}")
                    range_df = range_df[
                        (range_start-additional_expansion <= range_df[col_dict['position']]) & 
                        (range_df[col_dict['position']] <= range_end+additional_expansion)]

                    # file_name = f'{lead_SNP_id}-chr{chrom}.tsv.gz'
                    if len(range_df[range_df[col_dict['pvalue']] < p_threshold]) < 2:
                        logging.info(f"Less than 2 significant variants in {lead_SNP_id} loci")
                        continue

            else:
                logging.info(f"LD block file not found for {lead_SNP_id}")

                range_df = range_df[(
                        max(0, pos - additional_expansion) <= range_df[col_dict['position']]) & (
                                        range_df[col_dict['position']] <= pos + additional_expansion)]
            ############################# To speed up ###########################
            ## If the SNPs in range_df overlapped with cluster_pos_dict[chrom], regard this range_df
            if len(set(range_df[col_dict['position']]) & set(chrom_cluster_pos_list)) > 0:
                continue
            #####################################################################

            chrom_list.append(chrom)
            lead_SNP_list.append(lead_SNP_id)
            logging.info(f"chrom: {chrom}")
            logging.info(f"lead SNP id:{lead_SNP_id}")
            count += 1
            chrom_cluster_pos_list = chrom_cluster_pos_list + list(range_df[col_dict['position']])
            cluster_pos_dict[chrom] = list(set(chrom_cluster_pos_list))
            positions_list.append(range_df[col_dict['position']].tolist())
            del range_df

        # print(f"length of chrom_list: {len(chrom_list)}, {chrom_list}")
        # print(f"length of lead_SNP_list: {len(lead_SNP_list)}, {lead_SNP_list}")
        # print(f"length of positions_list: {len(positions_list)}, {positions_list}")
        cluster_summary_df = pd.DataFrame(
            {'chrom': chrom_list, 'lead_SNP': lead_SNP_list, 'positions': positions_list})
        cluster_summary_df['phenotype_id'] = phenotype_id

        if os.path.exists(outputpath) and os.path.getsize(outputpath) > 0:
            mode = 'a'
            header = False
        else:
            mode = 'w'
            header = True

        # torus input file include columns variant_id、loc、zscore
        cluster_summary_df.to_csv(outputpath, sep=const.output_spliter, 
                                  mode=mode, header=header, index=False)



    def preprocess_qtl_type(self, qtltype):
        # if qtltype == 'qtl1':
        self.preprocess_molqtl(self.qtl1_type, self.qtl1_grouped_dir, 
                                self.qtl1_file, self.qtl1_col_dict, 
                                self.qtl1_sep, self.qtl1_output_report,
                                self.qtl1_p_threshold, self.qtl1_preprocesed_dir,
                                self.qtl1_LD_window)

    def preprocess_molqtl(self, qtl_type, qtl_grouped_dir,qtl_file, qtl_col_dict, qtl_sep, qtl_output_report,
                          qtl_p_threshold, qtl_preprocesed_dir, qtl_LD_window):
        # group by trait
        logging.info(f'Start to split {qtl_type} file by gene')
        utils.delete_dir(qtl_grouped_dir)
        Path(qtl_grouped_dir).mkdir(exist_ok=True, parents=True)
        print(f"qtl_grouped_dir: {qtl_grouped_dir}")
        utils.split_file_by_col_name_qtl(qtl_grouped_dir, qtl_file,
                                     qtl_col_dict['chrom'], qtl_col_dict['phenotype_id'],
                                     readonly_cols=qtl_col_dict.values(), sep=qtl_sep)
        logging.info(f'finish split {qtl_type} file')
        # drop na or no significant snp gene file, sort by position
        total_gene_file_count = 0
        good_result_count = 0
        chrom_list = []
        pheno_file_list = []
        postions_list = []
        qtl_lead_snp_list = []
        logging.info(f'start to filter gene {qtl_type} file')

        if Path(qtl_LD_window).exists():
            qtl_LD_window_df = pd.read_csv(qtl_LD_window, sep=const.output_spliter)
            qtl_LD_window_df = qtl_LD_window_df.iloc[:-1,:]
        else:
            qtl_LD_window_df = None
        for chrom_dir in os.listdir(f'{qtl_grouped_dir}'):
        # for chrom_dir in ['15','6']:
            if not chrom_dir.startswith('.'):
                for qtl_file in os.listdir(f'{qtl_grouped_dir}/{chrom_dir}'):
                    # if eqtl_file.upper().startswith(const.gene_id_prefix):
                    total_gene_file_count += 1
                    pheno_id = utils.get_qtl_pheno_name(qtl_file)
                    print(f"qtl_pheno_id: {pheno_id}")
                    qtl_lead_snp, high_risk_snp_df = self.__prepare_qtl_data(
                        f'{qtl_grouped_dir}/{chrom_dir}/{qtl_file}', pheno_id, qtl_type, qtl_col_dict, 
                           qtl_p_threshold, qtl_preprocesed_dir, qtl_LD_window, qtl_LD_window_df)
                    if high_risk_snp_df is not None:
                        good_result_count += 1
                        chrom_list.append(chrom_dir)
                        pheno_file_list.append(qtl_file)
                        qtl_lead_snp_list.append(qtl_lead_snp)
                        postions_list.append(high_risk_snp_df[qtl_col_dict['position']].tolist())
                        logging.debug(f'good {qtl_type} pheno {qtl_file}')
        eqtl_filter_result = pd.DataFrame(
            {'chrom': chrom_list, 'pheno_file': pheno_file_list, 'lead_snp': qtl_lead_snp_list, 'positions': postions_list})
        eqtl_filter_result.to_csv(qtl_output_report, sep=const.output_spliter, index=False)
        logging.info(
            f'finish filter gene {qtl_type} file, {total_gene_file_count} gene files, good result {good_result_count}')
        return eqtl_filter_result

    def __prepare_qtl_data(self, qtl_trait_file_path, pheno_id, qtl_type, qtl_col_dict, 
                           qtl_p_threshold, qtl_preprocesed_dir, qtl_LD_window, qtl_LD_window_df):
        try:
            qtl_trait_df = pd.read_table(qtl_trait_file_path, sep=const.column_spliter, header=0,
                                          usecols=qtl_col_dict.values(),
                                          dtype={qtl_col_dict['chrom']: 'category',
                                                 qtl_col_dict['position']: 'Int64'})
        except Exception as e:
            logging.error(f'error to prepare {qtl_type} data: {qtl_trait_file_path}')
            logging.error(f'exception: {e}')

        if qtl_col_dict['position'] not in set(qtl_trait_df.columns):
            qtl_trait_df[qtl_col_dict['position']] = qtl_trait_df[qtl_col_dict['variant_id']].apply(lambda x: int(str(x).split('_')[1]))
        if qtl_col_dict['chrom'] not in set(qtl_trait_df.columns):
            qtl_trait_df[qtl_col_dict['chrom']] = qtl_trait_df[qtl_col_dict['variant_id']].apply(lambda x: str(x).split('_')[0].strip('chr'))
        # if qtl_col_dict['variant_id'] not in set(qtl_trait_df.columns):
        #     logging.info(f"variant id {qtl_col_dict['variant_id']} not matching {qtl_type} columns")
        #     qtl_trait_df[qtl_col_dict['variant_id']] = \
        #         'chr' + qtl_trait_df[qtl_col_dict['chrom']].astype(str) + '_' + \
        #         qtl_trait_df[qtl_col_dict['position']].astype(str)+ '_' + \
        #         qtl_trait_df[qtl_col_dict['ref']].astype(str)+ '_' + \
        #         qtl_trait_df[qtl_col_dict['alt']].astype(str)
                
        qtl_trait_df[Processor.VAR_ID_COL_NAME] = qtl_trait_df[qtl_col_dict['variant_id']].apply(lambda x: '_'.join(str(x).split('_')[:4]))
        utils.clean_data(qtl_trait_df, dup_consider_subset=Processor.VAR_ID_COL_NAME, keep_dup='first')

        utils.drop_indel_snp(qtl_trait_df, qtl_col_dict['alt'], qtl_col_dict['ref'])

        # print(f"qtl_trait_df pheno_id: {pheno_id} xxxxx {qtl_trait_df}")

        if qtl_trait_df.empty:
            del qtl_trait_df
            utils.delete_file_if_exists(qtl_trait_file_path)
            logging.info(f'No required data in {qtl_type} file {qtl_trait_file_path}')
            return None, None
        # else:
        print(f"{qtl_type}_trait_df_exist: {pheno_id}")
        pval_filter_qtl_df = \
            utils.filter_data_frame_by_p_value(qtl_trait_df, qtl_p_threshold,
                                                qtl_col_dict['pvalue'], inplace=False)
        pval_filter_qtl_df = pval_filter_qtl_df.sort_values(qtl_col_dict['pvalue'])
        # pval_filter_caqtl_df = qtl_trait_df.drop(
        #     qtl_trait_df[(qtl_trait_df[self.qtl_col_dict['pvalue']] == 0.0)].index,
        #     inplace=False)
        # logging.info(f"************qtl_p_threshold: {self.config_holder.qtl_p_threshold}")
        if pval_filter_qtl_df.empty:
            del qtl_trait_df
            del pval_filter_qtl_df
            utils.delete_file_if_exists(qtl_trait_file_path)
            return None, None
        # else:
        print(f"pval_filter_qtl_df_exist")

        # if qtl_LD_window_df is not None:
        #     print("qtl_LD_window_df exist")
        #     qtl_pheno_df = qtl_LD_window_df[qtl_LD_window_df['phenotype_id'] == pheno_id]
        #     chrom_list = qtl_pheno_df['chrom'].tolist()
        #     lead_SNP_list = qtl_pheno_df['lead_SNP'].tolist()
        #     positions_list = qtl_pheno_df['positions'].apply(lambda x: ast.literal_eval(x)).tolist()
        #     cluster_pos_dict = dict(zip(chrom_list, positions_list))
        #     count = len(lead_SNP_list)
        #     # phenotype_id_list = qtl_LD_window_df['phenotype_id'].tolist()
        # else:
        #     chrom_list = []
        #     lead_SNP_list = []
        #     positions_list = []
        #     cluster_pos_dict = {}
        #     count = 0

        # # if self.genomic_window == 'combined_LD_based_window':

        # self.__qtl_ld_interval(pval_filter_qtl_df, qtl_trait_df, 
        #                     self.window_size, 
        #                     qtl_preprocesed_dir,
        #                     qtl_col_dict, chrom_list, 
        #                     lead_SNP_list, positions_list, 
        #                     cluster_pos_dict, count, 
        #                     qtl_p_threshold,
        #                     qtl_LD_window, pheno_id,
        #                     additional_expansion = self.LD_additional_expansion
        #                     )
        # logging.info(f"qtl_trait_df: {qtl_trait_df}")
        qtl_lead_snp = list(qtl_trait_df[Processor.VAR_ID_COL_NAME])[0]
        qtl_trait_df.sort_values([qtl_col_dict['chrom'], qtl_col_dict['position']],
                                    inplace=True)
        qtl_trait_df.to_csv(qtl_trait_file_path, sep=const.output_spliter, index=False)
        logging.info(f"__prepare_qtl_data finishpheno: {pheno_id}")
        logging.info(f"qtl_lead_snp: {qtl_lead_snp}")
        return qtl_lead_snp, qtl_trait_df 



if __name__ == '__main__':
    processor = Processor()
    processor.preprocess_gwas()
    processor.preprocess_eqtl()
