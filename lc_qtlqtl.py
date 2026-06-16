import concurrent
import logging
import os
import sys
import traceback
import uuid
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from collections import defaultdict
from pathlib import Path
import pandas as pd
import numpy as np
import common.config
import common.constants
from xQTL_xQTL.api import run_tools_xqtl_xqtl_api
import json
from importlib import metadata
import pyranges as pr
from common import utils, global_data_process as gdp
from figures import report_data_processor as redp
from ranking import scoring as sc
from common import integrate_results

TOOL_SIG_COL_INFO = {"coloc": "overall_H4",
                     "ecaviar": "clpp",
                     "fastenloc": "GRCP",
                     "smr": "p_SMR"}


def __before_run_fastenloc_tools_check(global_config):
    logging.info(f'start check fastenloc')
    # utils.check_file_or_path_exist(global_config['input']['ld_block_loci_file'])
    # utils.check_file_or_path_exist(global_config['input']['eqtl_finemapping_file'])


def __before_run_coloc_tools_check(global_config):
    logging.info(f'start check coloc')


def __before_run_predixcan_tools_check(global_config):
    logging.info(f'start check predixcan')
    _model_db_path, _prediction_snp_covariance_path = utils.get_predixcan_ref_files(global_config)
    logging.info(f'Predixcan model db path: {_model_db_path}')
    logging.info(f'Predixcan SNP covariance path: {_prediction_snp_covariance_path}')
    # utils.check_file_or_path_exist(_model_db_path)
    # utils.check_file_or_path_exist(_prediction_snp_covariance_path)


def __before_run_smr_tools_check(global_config):
    logging.info(f'start check smr')
    utils.check_file_or_path_exist(global_config['input']['genecode'])
    __check_vcf(global_config)


def __before_run_ecaviar_tools_check(global_config):
    # no ref config to check
    logging.info(f'start check eCAVIAR')


def __before_run_twas_tools_check(global_config):
    logging.info(f'start check FUSION')
    twas_pos_path = utils.get_twas_ref_files(global_config)
    logging.info(f'TWAS pos file path: {twas_pos_path}')
    # utils.check_file_or_path_exist(twas_pos_path)
    __check_vcf(global_config)


# def __before_run_susie_tools_check(global_config):
#     logging.info(f'start check susie')


def __before_run_heels_tools_check(global_config):
    logging.info(f'start check heels')

def __before_run_hess_tools_check(global_config):
    logging.info(f'start check hess')

def __check_vcf(global_config):
    population = global_config.get('population', 'EUR').upper()
    ref_vcf_dir = global_config['input']['vcf']
    for chromosome in range(1, 23):
        input_vcf = os.path.join(ref_vcf_dir, population, f'chr{chromosome}.vcf.gz')
        utils.check_file_or_path_exist(input_vcf, False)


tools_func_map = {
    'coloc': {
        'check_fun': __before_run_coloc_tools_check,
        'run_fun_gwas_xqtl': run_tools_xqtl_xqtl_api.__preprocess_and_run_coloc,
    },
    }


################################################################################
#                                 Start analyses                               #
################################################################################

def run(config_file=None, log_file=None, parallel=False, tools_config=None, no_report=False, keep_intermediate=False):

    conda_env = os.getenv("CONDA_DEFAULT_ENV", "Unknown")
    logging.info(f"Current Conda environment: {conda_env}")
    cfg_list = []
    # retrieve config and study info
    if config_file is None:
        # default parameters
        study = common.constants.default_study
        cfg_list.append(common.constants.default_config)
    elif os.path.isdir(config_file):
        # specify study
        if config_file.endswith(os.sep):
            config_file = config_file.rstrip(os.sep)
        study = os.path.basename(config_file)
        for dir_path, _, file_names in os.walk(config_file):
            for cfg in file_names:
                if cfg.startswith('.'):
                    continue
                if not cfg.endswith('.yml') and not cfg.endswith('.yaml'):
                    continue
                cfg_list.append(os.path.join(dir_path, cfg))
    else:
        study = common.constants.default_study
        cfg_list.append(config_file)
    if len(cfg_list) == 0:
        logging.error(f'No config files to run.')
        return
        # loop to run all configs
    if log_file is None:
        log_file = f'{uuid.uuid4().hex}.log'

    report_list = []

    integrate_result_dict = __nested_dict()

    qtl_type_ls = []
    biological_context_ls = []
    total_numof_analyses = len(cfg_list)
    current_analysis_order = 0
    sub_study = 0
    for cfg in cfg_list:
        current_analysis_order = current_analysis_order + 1
        config_holder = common.config.ConfigHolder(single_config_file=cfg, 
                                                   study=study, parallel=parallel,
                                                   tools_config_file=tools_config,
                                                   keep_intermediate=keep_intermediate)
        Path(config_holder.study_dir).mkdir(parents=True, exist_ok=True)
        __init_logger(os.path.join(config_holder.study_dir, f'{log_file}'))

        ########################################################################
        #       1. Single study with specified QTL and Biological context      #
        ########################################################################

        population = config_holder.population.upper()
        biological1_context = config_holder.qtl1_biological_context
        qtl1_type = config_holder.qtl1_type.lower()

        # single_cfg_ensemble_result, integrate_result_sfg_output_file = \
        results_dict = \
            __run_single_cfg(config_holder, 
                             report_list, 
                             parallel, 
                             study, 
                             current_analysis_order, 
                             total_numof_analyses)
        sub_study = sub_study + 1
        calculated_schedule = int(80/total_numof_analyses * (current_analysis_order - 1))
        if os.path.exists('/process/'):
            with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
                schedule.write(str(calculated_schedule))
        else:
            pass
        schedule.close()


        integrate_result_dict[f"cfg_{sub_study}"] = {}

        try:
            utils.delete_dir(os.path.join(config_holder.qtl1_tool_parent_dir, 'vcf'), config_holder.keep_intermediate)
            utils.delete_dir(os.path.join(config_holder.qtl1_tool_parent_dir, 'vcf'), config_holder.keep_intermediate)

        except:
            logging.warning(f'failed to clean {config_holder.qtl1_tool_parent_dir}')

    # colocboost multitrait analyses gwas-qtls (GTEx)
    # tools_func_map['colocboost']['run_fun_gwas_xqtl'](
    #                     cfg_list, study, parallel, tools_config, keep_intermediate
    #                 )

        
    if os.path.exists('/process/'):
        with open(f"{os.path.join('/process/', 'process_schedule.log')}", 'w') as schedule:
            schedule.write(str(95))
    else:
        pass
    genecode_file = config_holder.global_config['input']['genecode']
    logging.info(f"genecode_file : {genecode_file}")
    lc3outpath = os.path.join(config_holder.final_json_dir,
                    f'locuscompare3_integrate_{datetime.now().strftime("%Y%m%d%H%M%S")}.json')
    logging.info(f"integrate_result_dict: {integrate_result_dict}")
    # results_collection(integrate_result_dict, lc3outpath, genecode_file)
    print(f"final_json_dir: {config_holder.final_json_dir}")
    it = integrate_results.INTEGRATE(genecode_file=genecode_file)  # 也可定制阈值/列名等
    it.run_all(integrate_result_dict,output_dir=config_holder.final_json_dir)

    if len(report_list) == 0:
        logging.warning(f'No results of specified tools found')
        return
    if no_report is False:
        redp.report_data_process(report_list)

    schedule.close()



def __nested_dict():
    return defaultdict(__nested_dict)

# def __caqtl_nearest_gene(row, genecode_df):

#     chr_match = str(row['chrom'])
#     genes_chr = genecode_df[genecode_df['chr'] == chr_match]
#     if genes_chr.empty:
#         return None
#     idx_min = (genes_chr['position'] - row['peak_center']).abs().idxmin()
#     return genes_chr.loc[idx_min, 'gene_id']



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
    genecode_df.loc[:, 'position'] = pd.NA
    genecode_df['position'].mask(genecode_df['strand'] == '+', genecode_df['start'], inplace=True)
    genecode_df['position'].mask(genecode_df['position'].isna(), genecode_df['end'], inplace=True)
    genecode_df.drop(columns=['start', 'end'], inplace=True)
    genecode_df['chr'] = genecode_df['chr'].apply(lambda x: x.strip('chr'))
    return genecode_df


def __init_logger(logfile):
    log_format = '%(asctime)s::%(levelname)s::%(name)s::%(filename)s::%(lineno)d::%(message)s'
    stdout_hd = logging.StreamHandler(sys.stdout)
    stdout_hd.setLevel(logging.INFO)
    logging.root.handlers = []
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.FileHandler(logfile),
            stdout_hd,
        ]
    )


def get_tools_path(dir_path, tool_name):
    for path_name in os.listdir(os.path.join(dir_path, tool_name, 'analyzed')):
        if os.path.isfile(os.path.join(dir_path, tool_name, 'analyzed', path_name)) and (
                path_name.startswith(tool_name) or path_name.startswith('report')):
            return os.path.join(dir_path, tool_name, 'analyzed', path_name)

def union_result_gene(tissues_ensemble_result_ls):
    print(f'union result gene: {tissues_ensemble_result_ls}')
    # 把每个tissue的基因统一
    all_tissue_df = pd.DataFrame()
    for single_cfg_ensemble_result in tissues_ensemble_result_ls:
        tmp = pd.read_csv(single_cfg_ensemble_result,sep='\t',usecols=['phenotype_id'])
        all_tissue_df = pd.concat([all_tissue_df,tmp],axis=0)
    all_tissue_df = all_tissue_df.dropna()
    all_gene = pd.DataFrame(list(set(all_tissue_df['phenotype_id'])))
    all_gene.columns = ['ensemble_gene_id']
    all_gene.index = all_gene['ensemble_gene_id']

    for single_cfg_ensemble_result in tissues_ensemble_result_ls:
        tmp = pd.read_csv(single_cfg_ensemble_result,sep='\t')
        tmp.index = tmp['phenotype_id']
        tmp = pd.concat([tmp,all_gene],axis=1)
        tmp['phenotype_id'] = tmp['ensemble_gene_id']
        tmp = tmp.iloc[:,:-1]
        tmp = tmp.reset_index(drop=True)
        tmp['geo_ranking'] = tmp.apply(lambda x: int(x.name) if pd.isna(x['geo_ranking']) else int(x['geo_ranking']), axis=1)
        tmp['avg_ranking'] = tmp.apply(lambda x: int(x.name) if pd.isna(x['avg_ranking']) else int(x['avg_ranking']), axis=1)
        tmp.to_csv(single_cfg_ensemble_result,sep='\t',index=False)
    pass



# def output_gwas_info(num_of_sig_gwas_SNP, gwas_cluster, extraFile_path):
#     gwas_cluster_df = pd.read_csv(gwas_cluster, sep='\t')
#     num_of_gwas_loci = len(gwas_cluster_df)
#     with open(f"{extraFile_path}", 'w') as output:
#         output.write(f"num_of_gwas_loci\t{num_of_gwas_loci}")
#     output.close()


# def __run_single_cfg(tools_param_list, config_holder, report_list, parallel, study):
def __run_single_cfg(config_holder, report_list, parallel, study, 
                     current_analysis_order, total_numof_analyses):
    start_time = datetime.now()
    tools_param_list = list(config_holder.global_config['tools'])
    qtl1_type = config_holder.qtl1_type.lower()
    qtl1_biological_context = config_holder.qtl1_biological_context
    
    logging.info(f'run tools_list: {tools_param_list}, '
                 f'qtl_type: {qtl1_type}, '
                 f'biological_context: {qtl1_biological_context}, '
                 f'start time: {start_time}')
    ############################################################################
    #                    1.1 Check QTL1 and QTL2 entry file exist               #
    ############################################################################
    utils.check_file_or_path_exist(config_holder.root_work_dir)

    utils.check_file_or_path_exist(config_holder.qtl1_file)
    utils.check_file_or_path_exist(config_holder.qtl2_file)

    # utils.check_file_or_path_exist(global_config['input']['vcf'])

    ############################################################################
    #                    1.2 Check tool require file exist                     #
    ############################################################################
    logging.info(type(tools_param_list))
    actually_tools_list = []
    if 'all' in tools_param_list:
        actually_tools_list = tools_func_map.keys()
    else:
        for tool in tools_param_list:
            logging.info(tool)
            if tools_func_map.get(tool):
                actually_tools_list.append(tool)
            else:
                logging.error(f'The {tool} tool is not recognized')
    logging.info("*****actually_tools_list*****")
    logging.info(actually_tools_list)

    smr_schedule = False
    ecaviar_schedule = False
    coloc_schedule = False
    fastenloc_schedule = False
    fusion_schedule = False
    predixcan_schedule = False

    for tool in actually_tools_list:
        if tool == 'smr':
            smr_schedule = True
            break
        if tool == 'ecaviar':
            ecaviar_schedule = True
            break
        if tool == 'fusion':
            fusion_schedule = True
            break
        if tool == 'coloc':
            coloc_schedule = True
            break
        if tool == 'fastenloc':
            fastenloc_schedule = True
            break
        if tool == 'predixcan':
            predixcan_schedule = True
            break

    config_schedule = {'smr': smr_schedule,
                        'ecaviar': ecaviar_schedule,
                        'coloc': coloc_schedule,
                        'fastenloc': fastenloc_schedule,
                        'fusion': fusion_schedule,
                        'predixcan': predixcan_schedule,
                        'heels': False,
                        # 'susie': False,
                        'hess': False}

    for tool in actually_tools_list:
        tools_func_map[tool]['check_fun'](config_holder.global_config)

    processor = gdp.Processor(config_holder)
    print(f"processor.rsidvcf: {processor.rsidvcf}")
    ############################################################################
    #                            1.3 Preprocess QTL                            #
    ############################################################################

    if not utils.file_exists(processor.qtl1_output_report):
        processor.preprocess_qtl_type('qtl1')

        qtl1_sig_df = pd.read_table(config_holder.qtl1_output_report, nrows=2)
        if qtl1_sig_df.shape[0] == 0:
            logging.warning(f'No significant records found in QTL file {config_holder.qtl1_file} '
                            f'by threshold {config_holder.qtl1_p_threshold}')
            return

    if not utils.file_exists(processor.qtl2_output_report):
        processor.preprocess_qtl_type('qtl2')

        qtl2_sig_df = pd.read_table(config_holder.qtl2_output_report, nrows=2)
        if qtl2_sig_df.shape[0] == 0:
            logging.warning(f'No significant records found in QTL file {config_holder.qtl2_file} '
                            f'by threshold {config_holder.qtl2_p_threshold}')
            return

    

    ############################################################################
    #                   1.5 Check preprocess file is exist                     #
    ############################################################################

    utils.check_path_exist_and_has_size(processor.qtl1_output_report)
    utils.check_path_exist_and_has_size(processor.qtl2_output_report)

    # ############################################################################
    # #                 1.6 Identify matching loci for GWAS/QTL                  #
    # ############################################################################

    # if set(['coloc','fastenloc','ecaviar','smr']) & set(actually_tools_list): 
    #         print(f"identifymatchinglociforgwasqtl")
    #         # processorloci = gwas_qtl_loci.ProcessGWASQTLLoci()  
    #         # processorloci.run(processor)


    results = dict()
    ############################################################################
    #                        1.7 Define rank output file                       #
    ############################################################################
    # rank_output_file = os.path.join(processor.qtl1_rank_dir,
    #                                 f'ensemble_ranking_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv')
    # integrate_result_sfg_output_file = os.path.join(processor.qtl1_rank_dir,
    #                                 f'ensemble_result_{datetime.now().strftime("%Y%m%d%H%M%S")}.tsv')
    if parallel:
        # 1) split susie out (run last)
        run_last_tool = "susie"
        parallel_tools = [t for t in actually_tools_list if t != run_last_tool]
        run_last = run_last_tool in actually_tools_list

        logging.info(f"Parallel running {parallel_tools}" + (f", then run {run_last_tool} last" if run_last else ""))

        results = {}
        exceptions = []

        # 2) parallel run (excluding susie)
        if parallel_tools:
            with ProcessPoolExecutor(max_workers=len(parallel_tools)) as executor:
                tool_futures = {}

                for tool in parallel_tools:
                    fun_key = "run_fun_gwas_xqtl"
                    fut = executor.submit(
                        tools_func_map[tool][fun_key],
                        processor,
                        current_analysis_order,
                        total_numof_analyses,
                        whether_schedule=config_schedule.get(tool, False),
                    )
                    tool_futures[fut] = tool

                for future in concurrent.futures.as_completed(tool_futures):
                    current_tool = tool_futures[future]
                    try:
                        results[current_tool] = future.result()
                        logging.info(f"Tool {current_tool} completed successfully!")
                    except Exception as exc:
                        exceptions.append((current_tool, exc))
                        logging.info(f"Tool {current_tool} failed with error: {exc}!")

        # 3) if any parallel tool failed, log + raise (before running susie)
        if exceptions:
            for tool, ex in exceptions:
                logging.error(f"[{tool}] " + "".join(traceback.TracebackException.from_exception(ex).format()))
            raise Exception([ex for _, ex in exceptions])

        # # 4) run susie last (serial)
        # if run_last:
        #     try:
        #         fun_key = "run_fun_gwas_xqtl"
        #         results[run_last_tool] = tools_func_map[run_last_tool][fun_key](
        #             processor,
        #             current_analysis_order,
        #             total_numof_analyses,
        #             whether_schedule=config_schedule.get(run_last_tool, False),
        #         )
        #         logging.info(f"Tool {run_last_tool} completed successfully (ran last)!")
        #     except Exception as exc:
        #         logging.error(f"Tool {run_last_tool} failed with error: {exc}!")
        #         logging.error("".join(traceback.TracebackException.from_exception(exc).format()))
        #         raise
    else:
        for tool in actually_tools_list:
            logging.info(f'Tool {tool} start running')
            try:
                report_file_path = tools_func_map[tool][f"run_fun_gwas_xqtl"](processor, current_analysis_order, total_numof_analyses, whether_schedule = config_schedule[tool])
                logging.info(f'Tool {tool} {report_file_path} completed successfully!')
            except Exception as error:
                logging.info(f'Tool {tool} failed with error: {error}!')
                raise error
            results[tool] = report_file_path
            # report_list.append({'trait': config_holder.gwas_trait, 'tool_name': tool, 'study1': study,
            #                     'biological_context': biological_context, 'report_path': report_file_path,
            #                     'cfg_pro': processor, 'rank_output_file': rank_output_file})
    # sc.run_ranking(rpt_obj=results, output_file_path=rank_output_file, out_integrate_results = integrate_result_sfg_output_file, 
    #                sample_size=processor.global_config['input']['gwas']['sample_size'], qtl_type=qtl_type)
    # logging.info(f'coloctools complete at: {datetime.now()},duration: {datetime.now() - start_time}')
    # logging.info(f'biological_context: {biological_context},trait: {config_holder.gwas_trait} coloctools complete')
    # return rank_output_file, integrate_result_sfg_output_file, qtl_type, biological_context

    if results is None or len(results) == 0:
        logging.warning(f'No results of specified tools found')
    else:
        pass

    return results


if __name__ == '__main__':
    start_time = datetime.now()
    logging.info(f'start run all coloctools, start time: {start_time}')
    parse_args = utils.parse_parameters()
    # run(parse_args.config_file, parse_args.tools_list, parse_args.log_file, parse_args.parallel,
    #     parse_args.tools_config, parse_args.no_report)
    # run(parse_args.config_file, parse_args.log_file, parse_args.parallel,

    run(parse_args.config_file, parse_args.log_file, True,
        parse_args.tools_config, parse_args.no_report, parse_args.keep_intermediate)
    logging.info(f'all coloctools complete at: {datetime.now()},duration: {datetime.now() - start_time}')


