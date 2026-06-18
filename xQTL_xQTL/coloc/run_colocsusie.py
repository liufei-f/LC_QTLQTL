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
from glob import glob
from importlib import metadata
import pyranges as pr

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
from common import utils, global_data_process as gdp, constants as const

def outputschedule(rownum, totalnum, current_analysis_order, total_numof_analyses, rank_dir):
    calculated_schedule = int(rownum/totalnum * 80/total_numof_analyses + 80/total_numof_analyses * (current_analysis_order - 1))
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    else:
        with open(f"{os.path.join(rank_dir, 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(calculated_schedule))
    schedule.close()

class COLOCSUSIE:
    COLOC_TOOL_NAME = 'susie'

    def __init__(self):
        logging.info(f'init qtl susie')

    def run(self, processor, input_files, 
            coloc_dir_input, coloc_dir_output,
            susie_input_dir, susie_output_dir, 
            final_report_dir,
            current_analysis_order = None, 
            total_numof_analyses = 1, whether_schedule = False):


        Path(susie_input_dir).mkdir(exist_ok=True, parents=True)
        Path(susie_output_dir).mkdir(exist_ok=True, parents=True)

        self.qtl1_type = processor.qtl1_type
        self.qtl2_type = processor.qtl2_type
        tools_config_file = processor.tools_config_file

        qtl1_sample_size = processor.qtl1_sample_size
        qtl2_sample_size = processor.qtl2_sample_size


        _qtl1data_type = processor._qtl1_type
        _qtl2data_type = processor._qtl2_type
        qtl1_threshold = processor.config_holder.qtl1_p_threshold
        qtl2_threshold = processor.config_holder.qtl2_p_threshold
        genecode_file = processor.global_config['input']['genecode']

        self.qtl1_col_dict = processor.qtl1_col_dict
        self.qtl2_col_dict = processor.qtl2_col_dict

        self.ref_vcf_dir = processor.ref_vcf_dir
        self.population = processor.population
        self.shell_command_plink_makebed_execute = 'plink --vcf {} --extract {} --make-bed --out {}'
        self.shell_command_plink_ld_matrix_execute = 'plink --bfile {} --r square --out {}'

        logging.info(f'init {self.qtl1_type} susie')
        logging.info(f'input_files: {input_files}')
        logging.info(f'coloc_dir_input: {coloc_dir_input}')
        logging.info(f'susie_input_dir: {susie_input_dir}')

        _p1, _p2, _p12 = self.__get_coloc_run_params(tools_config_file)
        # logging.info(f"_p1, _p2, _p12: {_p1}, {_p2}, {_p12}")

        logging.info(f'susie process_pheno')
        # genecode_df = self.__prepare_genecode(genecode_file)

        # gtex_config = processor.global_config['input'].get('gtex', {})
        # gtex_eqtl_expression_dir = gtex_config.get('gtex_eqtl_expression_dir')
        # gtex_eqtl_covariates_dir = gtex_config.get('gtex_eqtl_covariates_dir')
        # gtex_genotype_dir        = gtex_config.get('gtex_genotype_dir')

        # if all([gtex_eqtl_expression_dir, gtex_eqtl_covariates_dir, gtex_genotype_dir]):
        #     gene_info_file       = os.path.join(susie_input_dir, 'gene_info.tsv')
        #     susie_gtex_output_dir = os.path.join(susie_output_dir, 'gtex')
        #     self.run_susie_gtex(genecode_df, gtex_eqtl_expression_dir, gtex_eqtl_covariates_dir,
        #                         gtex_genotype_dir, gene_info_file, susie_gtex_output_dir)

        results = self.run_susie_qtl_qtl(input_files, coloc_dir_input, coloc_dir_output,
                                          susie_input_dir, susie_output_dir, 
                                          qtl1_sample_size, qtl2_sample_size,
                                          _qtl1data_type, _qtl2data_type, _p1, _p2, _p12,
                                          final_report_dir)
        utils.delete_dir(coloc_dir_output, processor.keep_intermediate)

        return results


    def run_susie_qtl_qtl(self, input_files, coloc_dir_input, coloc_dir_output,
                                          susie_input_dir, susie_output_dir, 
                                          qtl1_sample_size, qtl2_sample_size,
                                          _qtl1data_type, _qtl2data_type, _p1, _p2, _p12,
                                          final_report_dir):
        logging.info(f"running susie {coloc_dir_input}")
        start_time = datetime.now()
        for each in input_files:
            
            qtl1_coloc_path = each
            qtl1_coloc_filename = os.path.basename(each)
            
            qtl2_coloc_filename = 'qtl2_'+'_'.join(qtl1_coloc_filename.split('_')[1:])
            qtl2_coloc_path = os.path.join(coloc_dir_input, qtl2_coloc_filename)

            qtl1_colocsusie_path = os.path.join(susie_input_dir, qtl1_coloc_filename)
            qtl2_colocsusie_path = os.path.join(susie_input_dir, qtl2_coloc_filename)

            chrom = qtl1_coloc_filename.split('_')[1].strip('chr')
            qtl1 = pd.read_csv(qtl1_coloc_path, sep='\t')
            qtl1.index = qtl1['var_id_']
            qtl2 = pd.read_csv(qtl2_coloc_path, sep='\t')
            qtl2.index = qtl2['var_id_']


            logging.info(f"qtl1_coloc_path: {qtl1_coloc_path}")
            logging.info(f"qtl2_coloc_path: {qtl2_coloc_path}")
            # min_pos = str(min(gwas[self.gwas_col_dict['position']]))
            # max_pos = str(max(gwas[self.gwas_col_dict['position']]))


            qtl1_phenotype_id = (
                qtl1_coloc_filename
                .removesuffix('.tsv.gz')
                .split(f'_{self.qtl2_type}_')[0]
                .removeprefix(f"qtl1_chr{chrom}_{self.qtl1_type}_")
            )
            qtl2_phenotype_id = (
                qtl2_coloc_filename
                .removesuffix('.tsv.gz')
                .split(f'_{self.qtl2_type}_')[1]
            )
            logging.info(f"check phenotype IDs: {qtl1_phenotype_id}, {qtl2_phenotype_id}")

            filename = f'chr{chrom}_{self.qtl1_type}_{qtl1_phenotype_id}_{self.qtl2_type}_{qtl2_phenotype_id}'
            plinkld_outname = os.path.join(susie_input_dir, filename)
            output_file = os.path.join(susie_output_dir, f'{filename}.tsv')

            rscript_path = os.path.join(os.path.dirname(Path(__file__).resolve()), 'rscript', 'colocsusie_standard.R')
            # logging.info(f"{rscript_path}")
            common_variants = set(qtl1.index).intersection(set(qtl2.index))
            pd.DataFrame(list(common_variants)).to_csv(os.path.join(susie_input_dir, f'{filename}.snplist.txt'),index=False,header=False)
            os.system(self.shell_command_plink_makebed_execute.format(os.path.join(self.ref_vcf_dir, self.population, f"chr{chrom}.vcf.gz"),
                                                                        os.path.join(susie_input_dir, f'{filename}.snplist.txt'),
                                                                        plinkld_outname
                                                                        ))
            os.system(self.shell_command_plink_ld_matrix_execute.format(plinkld_outname,
                                                                        plinkld_outname
                                                                        ))
            LD_bim = pd.read_csv(f'{plinkld_outname}.bim',sep='\t',header=None)
            LD_matrix = pd.read_csv(f'{plinkld_outname}.ld',sep='\t',header=None)
            
            logging.info(f"susie_plinkld_outname: {plinkld_outname}")
            LD_matrix.index = LD_bim[1]
            # LD_matrix = LD_matrix.dropna(axis=0, how='all')
            # LD_matrix.columns = LD_bim[1]

            qtl1 = qtl1.loc[list(LD_matrix.index)]
            qtl2 = qtl2.loc[list(LD_matrix.index)]
            # LD_matrix = LD_matrix.loc[common_variants]

            qtl1['z'] = qtl1['beta'] / qtl1['se']
            qtl2['z'] = qtl2['beta'] / qtl2['se']
            qtl1['n'] = qtl1_sample_size
            qtl2['n'] = qtl2_sample_size
            qtl1 = qtl1[[ 'z', 'n' , 'var_id_']]
            qtl1.columns = ['z', 'n', 'variant']
            qtl2 = qtl2[[ 'z', 'n' , 'var_id_']]
            qtl2.columns = ['z', 'n', 'variant']
            qtl2.to_csv(qtl2_colocsusie_path,sep='\t',index=False)
            qtl1.to_csv(qtl1_colocsusie_path,sep='\t',index=False)
            LD_matrix.to_csv(f'{plinkld_outname}.ld',sep='\t',header=None)

            logging.info(f'Rscript --no-save --no-restore {rscript_path} '
                    f'{output_file} {qtl1_colocsusie_path} {qtl2_colocsusie_path} '
                    f'{qtl1_sample_size} {qtl2_sample_size} {_qtl1data_type} {_qtl2data_type} '
                    f'{qtl1_phenotype_id} {qtl2_phenotype_id} '
                    f'{_p1} {_p2} {_p12} {plinkld_outname}.ld')
            
            os.system(f'Rscript --no-save --no-restore {rscript_path} '
                    f'{output_file} {qtl1_colocsusie_path} {qtl2_colocsusie_path} '
                    f'{qtl1_sample_size} {qtl2_sample_size} {_qtl1data_type} {_qtl2data_type} '
                    f'{qtl1_phenotype_id} {qtl2_phenotype_id} '
                    f'{_p1} {_p2} {_p12} {plinkld_outname}.ld')

            # if Path(output_file).exists():
            #     logging.info("{} completed".format(output_file))
            #     coloc_out_file = os.path.join(coloc_dir_output, f'{filename}.tsv')
            #     self.integrate_coloc_susie(coloc_out_file, output_file)
            # else:
            #     logging.info("SUSIE {} FAILED".format(output_file))

        final_output_file = self.get_output_file(final_report_dir)
        self.__analyze_result(susie_output_dir, final_output_file)

        if not os.path.exists(final_output_file) or os.path.getsize(final_output_file) <= 0:
            logging.warning(f'Process completed, duration {datetime.now() - start_time}, no result found')
        else:
            logging.info(
                f'Process completed, duration {datetime.now() - start_time}, with params p1: {_p1} p2:{_p2} p12:{_p12}, check {final_output_file} for result!')

        return final_output_file


    def __get_coloc_run_params(self, tools_config):
        params = utils.get_tools_params_dict(self.COLOC_TOOL_NAME, tools_config)
        _p1 = 1.0E-4 if params.get('p1') is None else params['p1']
        _p2 = 1.0E-4 if params.get('p2') is None else params['p2']
        _p12 = 1.0E-5 if params.get('p12') is None else params['p12']
        return _p1, _p2, _p12

    def integrate_coloc_susie(self, coloc_out_file, susie_out_file):
        coloc_df = pd.read_csv(coloc_out_file, sep='\t')
        susie_df = pd.read_csv(susie_out_file, sep='\t')
        coloc_df.index = coloc_df['snp']
        coloc_susie_df = pd.merge(
                            coloc_df[['chrom','phenotype_id','lead_snp','SNP.PP.H4','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','overall_H4']],
                            susie_df, 
                            left_index=True,
                            right_index=True,
                            how='left')
        coloc_susie_df.rename(columns={'result.pip': 'susie_pip'}, inplace=True)
        coloc_susie_df.to_csv(susie_out_file, sep=const.output_spliter, header=True, index=False)

    def __analyze_result(self, output_dir, final_result_file):
        logging.info(f"output_dir: {output_dir}")
        logging.info(f"final_result_file: {final_result_file}")
        single_result_list = []
        for single_result in os.listdir(output_dir):
            if not (single_result.endswith('.tsv') or single_result.endswith('.tsv.gz')):
                continue
            df = pd.read_table(os.path.join(output_dir, single_result))
            single_result = single_result.replace('gwas_','').replace('.tsv.gz','').replace('.tsv','')
            phenotype_id = self.__get_phenotype_id_from_filename_diffqtl(single_result)
            
            df['phenotype_id'] = phenotype_id
            # filtered_df = df[df['hit1'] == df['hit2']]
            # filtered_df = filtered_df[['hit2', 'PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf', 'idx1', 'idx2', 'lead_snp', 'locus_range']]
            # filtered_df.rename(columns={'hit2': 'snp'}, inplace=True)

            single_result_list.append(df)
        if len(single_result_list) == 0:
            return
        report_df = pd.concat(single_result_list)
        # report_df.sort_values(by=['overall_H4', 'SNP.PP.H4'], ascending=False, inplace=True)
        # # report_df.rename(columns={'phenotype_id': 'gene_id'}, inplace=True)
        # report_df.drop_duplicates(subset=['snp', 'SNP.PP.H4', 'phenotype_id'], inplace=True)
        # report_df = report_df.round(4)
        report_df.to_csv(final_result_file, sep=const.output_spliter, header=True, index=False)

    def __get_phenotype_id_from_filename_diffqtl(self, filename):
        if self.qtl1_type == 'caqtl':
            return '_'.join(filename.split('_')[:3])
        else:
            return filename.split('_')[0]

    def get_output_file(self, working_dir):
        final_out_dir = os.path.join(working_dir, "analyzed")
        Path(final_out_dir).mkdir(parents=True, exist_ok=True)
        _output_file_name = f'{self.COLOC_TOOL_NAME}_output_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv.gz'
        return os.path.join(final_out_dir, _output_file_name)



    def __prepare_genecode(genecode_file):
        genecode_df = pr.read_gtf(genecode_file, as_df=True)
        genecode_df.drop(labels=genecode_df[genecode_df['Feature'] != 'gene'].index, inplace=True)
        genecode_df.reset_index(drop=True, inplace=True)
        genecode_df.drop(columns=[col for col in genecode_df if
                                    col not in ['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_name']],
                            inplace=True)
        genecode_df.rename({'Chromosome': 'chr', 'Start': 'start', 'End': 'end', 'Strand': 'strand'},
                            axis='columns', inplace=True)
        # gene_id_df = genecode_df['pheno_id'].str.split('.', n=1, expand=True)
        # genecode_df.loc[:, 'pheno_id'] = gene_id_df[0]
        # genecode_df.loc[:, 'position'] = pd.NA
        # genecode_df['position'].mask(genecode_df['strand'] == '+', genecode_df['start'], inplace=True)
        # genecode_df['position'].mask(genecode_df['position'].isna(), genecode_df['end'], inplace=True)
        # genecode_df.drop(columns=['start', 'end'], inplace=True)
        genecode_df['chr'] = genecode_df['chr'].apply(lambda x: x.strip('chr'))
        return genecode_df


if __name__ == '__main__':
    processor = gdp.Processor()
    if len(os.listdir(processor.gwas_cluster_output_dir)) == 0 or not os.path.exists(
            processor.qtl_output_report):
        raise ValueError(f'Dependant files not found, did you run global_data_process?')
    _working_dir = os.path.join(processor.tool_parent_dir, SUSIE.COLOC_TOOL_NAME)
    Path(_working_dir).mkdir(exist_ok=True, parents=True)
    _gwas_sample_size = processor.global_config['input']['gwas']['sample_size']
    _qtl_sample_size = processor.global_config['input']['qtl1'].get('sample_size', 948)
    _gwas_type = processor.global_config['input']['gwas'].get('type', 'cc')
    _qtl_type = processor.global_config['input']['qtl1'].get('type', 'quant')
    coloc = SUSIE()
    coloc.run(_working_dir,
              gdp.Processor.VAR_ID_COL_NAME,
              processor.gwas_cluster_output_dir,
              processor.gwas_col_dict,
              _gwas_sample_size,
              processor.qtl_output_report,
              processor.qtl_grouped_dir,
              processor.qtl_col_dict,
              _qtl_sample_size,
              _gwas_type,
              _qtl_type)
