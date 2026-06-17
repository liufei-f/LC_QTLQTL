import asyncio
import os
from pathlib import Path
import sys
import pandas as pd
from datetime import datetime

sys.path.append(
    os.path.abspath(os.path.join(os.path.join(os.path.dirname(Path(__file__).resolve()), os.pardir), os.pardir)))
print(sys.path)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print(sys.path)

import xQTL_xQTL.coloc.run_colocsusie as rcs
# from caQTL.predixcan import run_predixcan as rp, gwas_data_processor as pgdp
# from caQTL.twas import run_twas as rt
from common import global_data_process as gdp
from common import utils, constants as const
import logging
from glob import glob

import xQTL_xQTL.coloc.coloc_dataprocessor as coloc_dp

# import xQTL_xQTL.colocboost.dataprocessor as colocboost_dp


thresholds = {'coloc': 0.75,
             'ecaviar': 0.01,
             'fastenloc': 0.5}

metrics = {'coloc': 'overall_H4',
           'ecaviar': 'clpp',
           'fastenloc': 'GRCP',
           'smr': 'p_SMR'}

def sort_results(tool):
    if tool in ['coloc','ecaviar','fastenloc']:
        return False
    else:
        return True




######################
###    1. coloc    ###
######################

def __preprocess_and_run_coloc(processor, current_analysis_order, total_numof_analyses, whether_schedule):
    print(f"__preprocess_and_run_coloc: {processor.global_config['input']['qtl1']['qtl_type']}")

    tool = 'coloc'
    results = processor.results
    result_path = results.get(tool) if isinstance(results, dict) else None

    if not result_path or not Path(result_path).exists():
        logging.info(f"coloc")
        processorloci = coloc_dp.COLOC2QTLLOCI()  
        processorloci.run(processor)

        # coloc = rc.Coloc()
        if 'gtex' not in processor.global_config or 'LDsketch_meta_path' not in processor.global_config:
            colocsusie = rcs.COLOCSUSIE()
            coloc_dir_input = os.path.join(processor.coloc_base_dir, 'input')
            coloc_dir_output = os.path.join(processor.coloc_base_dir, 'output')
            input_files = glob(os.path.join(coloc_dir_input, 'qtl1*'))
            final_report_dir = processor.coloc_base_dir
            susie_input_dir = os.path.join(processor.susie_base_dir, 'input')
            susie_output_dir = os.path.join(processor.susie_base_dir, 'output')
            results = colocsusie.run(processor, input_files, 
                                                    coloc_dir_input, coloc_dir_output,
                                                    susie_input_dir, susie_output_dir, 
                                                    final_report_dir,
                                                    current_analysis_order = current_analysis_order, 
                                                    total_numof_analyses = total_numof_analyses, 
                                                    whether_schedule = whether_schedule)

        logging.info(f"coloc done")

        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        # validate_two_windowsizes(tool, combined_ld_window_results, fixed_window_results, final_report)
    else:
        final_report = processor.results[tool]

    return final_report


######################
###  2. fastenloc  ###
######################

def __preprocess_and_run_fastenloc(processor, current_analysis_order, total_numof_analyses, whether_schedule):
    print(f"__preprocess_and_run_fastenloc: {processor.global_config['input']['qtl1']['qtl_type']}")

    tool = 'fastenloc'
    results = processor.results
    result_path = results.get(tool) if isinstance(results, dict) else None

    if not result_path or not Path(result_path).exists():
        logging.info(f"{tool}")
        processorloci = fastenloc_dp.FASTENLOCGWASQTLLoci()  
        processorloci.run(processor)
        rf_obj = rf.Fastenloc()
        output_analyze_output_dir = os.path.join(processor.qtl1_fastenloc_base_dir, 'analyzed')
        output_fastenlocpreprocess_dir = os.path.join(processor.qtl1_fastenloc_base_dir, 'fastenloc_input')
        results = rf_obj.run(processor, output_analyze_output_dir, output_fastenlocpreprocess_dir)
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    else:
        final_report = processor.results[tool]

    return final_report

######################
###   3. ecaviar   ###
######################

def __preprocess_and_run_ecaviar(processor, current_analysis_order, total_numof_analyses, whether_schedule):
    print(f"__preprocess_and_run_ecaviar: {processor.global_config['input']['qtl1']['qtl_type']}")

    tool = 'ecaviar'
    results = processor.results
    result_path = results.get(tool) if isinstance(results, dict) else None

    if not result_path or not Path(result_path).exists():
        logging.info(f"{tool}")
        processorloci = ecaviar_dp.ECAVIARGWASQTLLoci()  
        processorloci.run(processor)
        ecaviar = run_e.ECaviar()
        candidate_data_dir = os.path.join(processor.qtl1_ecaviar_base_dir, 'candidate')
        outdir = os.path.join(processor.qtl1_ecaviar_base_dir, "analyzed")
        final_report_dir = processor.qtl1_ecaviar_base_dir
        results  = ecaviar.run(
            processor, candidate_data_dir, outdir, final_report_dir,
            current_analysis_order=current_analysis_order,
            total_numof_analyses=total_numof_analyses,
            whether_schedule=whether_schedule)
    else:
        final_report = processor.results[tool]
        
    logging.info(f"ecaviarresults: {final_report}")
    return final_report


######################
###     4. smr     ###
######################

def __preprocess_and_run_smr(processor, current_analysis_order, total_numof_analyses, whether_schedule):
    print(f"__preprocess_and_run_smr: {processor.global_config['input']['qtl1']['qtl_type']}")

    tool = 'smr'
    results = processor.results
    result_path = results.get(tool) if isinstance(results, dict) else None

    if not result_path or not Path(result_path).exists():
        logging.info(f"{tool}")
        # prepare ldref file finished
        if len(os.listdir(processor.gwas_output_dir)) == 0:
            logging.warning('Dependent files not found, did you run gwas_data_processor?')
            return None
        processorloci = smr_dp.SMRGWASQTLLoci()  
        processorloci.run(processor)
        smr = rs.Smr()
        coloc_dir_input = os.path.join(processor.qtl1_coloc_base_dir, 'input')
        input_dir = os.path.join(processor.qtl1_smr_base_dir, 'input')
        output_dir = os.path.join(processor.qtl1_smr_base_dir, 'output')
        smr_dir_ldref = os.path.join(processor.qtl1_smr_base_dir, 'ldref')
        vcf_output_dir = os.path.join(processor.qtl1_smr_base_dir, 'vcf')
        final_report_dir = processor.qtl1_smr_base_dir
        results =  smr.run(processor, input_dir, 
                                                output_dir, coloc_dir_input, 
                                                vcf_output_dir, smr_dir_ldref, 
                                                final_report_dir)

        logging.info(f"{tool} fixed_window")
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")


    else:
        final_report = processor.results[tool]
    
    return final_report


###    8. susie    ###
######################

def __preprocess_and_run_susie(processor, current_analysis_order, total_numof_analyses, whether_schedule):
    print(f"__preprocess_and_run_susie: {processor.global_config['input']['qtl1']['qtl_type']}")

    tool = 'susie'
    results = processor.results
    result_path = results.get(tool) if isinstance(results, dict) else None
    
    # if not result_path or not Path(result_path).exists():
    #     if processor.genomic_window == 'combined_LD_based_window' and processor.whether_fixed_window == True:
    #         logging.info(f"susie combined_LD_based_window")
    #         susie = rsu.SUSIE()
    #         coloc_dir_input = os.path.join(processor.qtl1_coloc_base_dir, 'input')
    #         coloc_dir_output = os.path.join(processor.qtl1_coloc_base_dir, 'output')
    #         susie_input_dir = os.path.join(processor.qtl1_combined_ld_susie_base_dir, 'input')
    #         susie_output_dir = os.path.join(processor.qtl1_combined_ld_susie_base_dir, 'output')
    #         input_files = glob(os.path.join(coloc_dir_input, 'gwas*'))
    #         final_report_dir = processor.qtl1_combined_ld_susie_base_dir
    #         results = susie.run(processor, input_files, coloc_dir_input, coloc_dir_output,
    #                                                susie_input_dir, susie_output_dir, final_report_dir,
    #                                                current_analysis_order = current_analysis_order, 
    #                                                total_numof_analyses = total_numof_analyses, 
    #                                                whether_schedule = whether_schedule)

    #         logging.info(f"susie fixed_window")
    #         susie = rsu.SUSIE()
    #         coloc_dir_input = os.path.join(processor.qtl1_fixed_win_coloc_base_dir, 'input')
    #         coloc_dir_output = os.path.join(processor.qtl1_fixed_win_coloc_base_dir, 'output')
    #         susie_input_dir = os.path.join(processor.qtl1_fixed_win_susie_base_dir, 'input')
    #         susie_output_dir = os.path.join(processor.qtl1_fixed_win_susie_base_dir, 'output')
    #         input_files = glob(os.path.join(coloc_dir_input, 'gwas*'))
    #         final_report_dir = processor.qtl1_fixed_win_susie_base_dir
    #         fixed_window_results = susie.run(processor, input_files, coloc_dir_input, coloc_dir_output,
    #                                         susie_input_dir, susie_output_dir, final_report_dir,
    #                                         current_analysis_order = current_analysis_order, 
    #                                         total_numof_analyses = total_numof_analyses, 
    #                                         whether_schedule = whether_schedule)
    #         timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    #         final_report = os.path.join(
    #             processor.qtl1_tool_parent_dir,
    #             f"{tool}_output_{timestamp}.tsv.gz"
    #         )
    #         print(f"Finalreport{tool} {final_report}")
    #         validate_two_windowsizes(tool, results, fixed_window_results, final_report)
    #     elif processor.genomic_window == 'fixed_GWAS_Loci_window':
    #         logging.info(f"susie fixed_window only")
    #         susie = rsu.SUSIE()
    #         coloc_dir_input = os.path.join(processor.qtl1_fixed_win_coloc_base_dir, 'input')
    #         coloc_dir_output = os.path.join(processor.qtl1_fixed_win_coloc_base_dir, 'output')
    #         susie_input_dir = os.path.join(processor.qtl1_fixed_win_susie_base_dir, 'input')
    #         susie_output_dir = os.path.join(processor.qtl1_fixed_win_susie_base_dir, 'output')
    #         input_files = glob(os.path.join(coloc_dir_input, 'gwas*'))
    #         final_report_dir = processor.qtl1_fixed_win_susie_base_dir
    #         final_report = susie.run(processor, input_files, coloc_dir_input, coloc_dir_output,
    #                                 susie_input_dir, susie_output_dir, final_report_dir,
    #                                 current_analysis_order = current_analysis_order, 
    #                                 total_numof_analyses = total_numof_analyses, 
    #                                 whether_schedule = whether_schedule)
    #     elif processor.genomic_window == 'combined_LD_based_window' and processor.whether_fixed_window == False:
    #         logging.info(f"susie combined_LD_based_window only")
    #         susie = rsu.SUSIE()
    #         coloc_dir_input = os.path.join(processor.qtl1_coloc_base_dir, 'input')
    #         coloc_dir_output = os.path.join(processor.qtl1_coloc_base_dir, 'output')
    #         susie_input_dir = os.path.join(processor.qtl1_combined_ld_susie_base_dir, 'input')
    #         susie_output_dir = os.path.join(processor.qtl1_combined_ld_susie_base_dir, 'output')
    #         input_files = glob(os.path.join(coloc_dir_input, 'gwas*'))
    #         final_report_dir = processor.qtl1_combined_ld_susie_base_dir
    #         final_report = susie.run(processor, input_files, coloc_dir_input, coloc_dir_output,
    #                                 susie_input_dir, susie_output_dir, final_report_dir,
    #                                 current_analysis_order = current_analysis_order, 
    #                                 total_numof_analyses = total_numof_analyses, 
    #                                 whether_schedule = whether_schedule)
    # else:
    #     final_report = processor.results[tool]

    # return final_report


######################
###    9. heels    ###
######################

def __preprocess_and_run_heels(processor, current_analysis_order, total_numof_analyses, whether_schedule):
    print(f"__preprocess_and_run_heels: {processor.global_config['input']['qtl1']['qtl_type']}")

    tool = 'heels'
    results = processor.results
    result_path = results.get(tool) if isinstance(results, dict) else None

    if not result_path or not Path(result_path).exists():
        if processor.genomic_window == 'combined_LD_based_window' and processor.whether_fixed_window == True:
            logging.info(f"{tool} combined_LD_based_window")
            heels = rheels.Heels()
            gwas_cluster_output_dir = processor.gwas_ldbased_cluster_output_dir
            gwas_cluster_summary = processor.gwas_ldbased_cluster_summary
            heels.run(processor, gwas_cluster_output_dir, gwas_cluster_summary)

            logging.info(f"{tool} fixed_window")
            heels = rheels.Heels()
            gwas_cluster_output_dir = processor.gwas_fixed_cluster_output_dir
            gwas_cluster_summary = processor.gwas_fixed_cluster_summary
            heels.run(processor, gwas_cluster_output_dir, gwas_cluster_summary)


        elif processor.genomic_window == 'fixed_GWAS_Loci_window':
            logging.info(f"{tool} fixed_window only")
            heels = rheels.Heels()
            gwas_cluster_output_dir = processor.gwas_fixed_cluster_output_dir
            gwas_cluster_summary = processor.gwas_fixed_cluster_summary
            heels.run(processor, gwas_cluster_output_dir, gwas_cluster_summary)

        elif processor.genomic_window == 'combined_LD_based_window' and processor.whether_fixed_window == False:
            logging.info(f"{tool} combined_LD_based_window only")
            heels = rheels.Heels()
            gwas_cluster_output_dir = processor.gwas_ldbased_cluster_output_dir
            gwas_cluster_summary = processor.gwas_ldbased_cluster_summary
            heels.run(processor, gwas_cluster_output_dir, gwas_cluster_summary)




######################
###    10. hess    ###
######################

def __preprocess_and_run_hess(processor, current_analysis_order, total_numof_analyses, whether_schedule):
    print(f"__preprocess_and_run_hess: {processor.global_config['input']['qtl1']['qtl_type']}")

    tool = 'hess'
    results = processor.results
    result_path = results.get(tool) if isinstance(results, dict) else None

    if not result_path or not Path(result_path).exists():
        if processor.genomic_window == 'combined_LD_based_window' and processor.whether_fixed_window == True:
            logging.info(f"{tool} combined_LD_based_window")
            hess = rhess.HESSRunner()
            gwas_cluster_summary = processor.gwas_ldbased_cluster_summary
            final_report = hess.run(processor, gwas_cluster_summary)

            logging.info(f"{tool} fixed_window")
            hess = rhess.HESSRunner()
            gwas_cluster_summary = processor.gwas_fixed_cluster_summary
            final_report = hess.run(processor, gwas_cluster_summary)


        elif processor.genomic_window == 'fixed_GWAS_Loci_window':
            logging.info(f"{tool} fixed_window only")
            hess = rhess.HESSRunner()
            gwas_cluster_summary = processor.gwas_fixed_cluster_summary
            final_report = hess.run(processor, gwas_cluster_summary)

        elif processor.genomic_window == 'combined_LD_based_window' and processor.whether_fixed_window == False:
            logging.info(f"{tool} combined_LD_based_window only")
            hess = rhess.HESSRunner()
            gwas_cluster_summary = processor.gwas_ldbased_cluster_summary
            final_report = hess.run(processor, gwas_cluster_summary)

    else:
        final_report = processor.results[tool]

    return final_report




if __name__ == '__main__':


    tool = 'coloc'
    fixed_window_results = '/home/users/nus/e1124850/scratch/lc3testgwas/processed/CONF/lc3test_GCST004610/plasma/EUR/pqtl/fixed_win/coloc/analyzed/coloc_output_20250812154323.tsv.gz'
    results = '/home/users/nus/e1124850/scratch/lc3testgwas/processed/CONF/lc3test_GCST004610/plasma/EUR/pqtl/combined_ld/coloc/analyzed/coloc_output_20250812155311.tsv.gz'
    final_report = '/home/users/nus/e1124850/scratch/lc3testgwas/processed/CONF/lc3test_GCST004610/plasma/EUR/pqtl/coloc_output_20250812155325.tsv.gz'
    ff= '/home/users/nus/e1124850/scratch/lc3testgwas/processed/CONF/lc3test_GCST004610/plasma/EUR/pqtl/coloc_output_20250812155325.tsv.gz'
    validate_two_windowsizes(tool, fixed_window_results, results, final_report)


