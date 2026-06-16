from pathlib import Path
import os
import sys
import common.constants
from common import utils
import logging


class ConfigHolder:

    def __init__(self, single_config_file=common.constants.default_config, study=common.constants.default_study,
                 parallel=False, tools_config_file=None, keep_intermediate=False):
        

        ########################################################################
        ##                          global setting                            ##
        ########################################################################

        if study is None:
            study = common.constants.default_study
        if single_config_file is None:
            single_config_file = common.constants.default_config

        self.global_config = utils.read_config(single_config_file)
        self.root_work_dir = self.global_config['working_dir']
        self.tools = self.global_config['tools']
        print(f"global_config: {self.global_config}")
        self.output_preprocessed_dir = os.path.join(self.root_work_dir, 'preprocessed')
        self.output_processed_dir = os.path.join(self.root_work_dir, 'processed')
        self.min_matching_number = self.global_config['min_matching_number']
        # Path(self.output_preprocessed_dir).mkdir(exist_ok=True, parents=True)
        # Path(self.output_processed_dir).mkdir(exist_ok=True, parents=True)     
        self.study_dir = os.path.join(self.output_processed_dir, study)

        self.genomic_window = self.global_config['input']['genomic_window']
        self.window_size = self.global_config['input']['window_size']
        self.LD_r2_filter = self.global_config['input']['LD_r2_filter']
        self.LD_additional_expansion = self.global_config['input']['LD_additional_expansion']        
        self.population = self.global_config['population']
        self.genecode_file = self.global_config['input']['genecode']

        self.report_path = os.path.join(self.study_dir) # trait/processed/default/

        self.ref_vcf_dir = self.global_config['input']['vcf']
        print(f"self.ref_vcf_dir: {self.ref_vcf_dir}")


        if 'ld_block_vcf_dir' in self.global_config['input'] and self.global_config['input']['ld_block_vcf_dir']:
            self.ld_block_vcf_dir = self.global_config['input']['ld_block_vcf_dir']
        else:
            self.ld_block_vcf_dir = False
        print(f"self.ld_block_vcf_dir: {self.ld_block_vcf_dir}")

        if ('vcf_rsid' in self.global_config['input']) and \
            self.global_config['input']['vcf_rsid']:
            self.rsidvcf = self.global_config['input']['vcf_rsid']
            print(f"self.global_config['input']: {self.global_config['input']}")
            print(f"self.rsidvcf: {self.rsidvcf}")
        else:
            self.rsidvcf = self.ref_vcf_dir
            print(f"no rsidvcf")
        print(f"self.rsidvcf: {self.rsidvcf}")

        self.parallel = parallel


        # tools parameter config file
        self.tools_config_file = tools_config_file
        self.keep_intermediate = keep_intermediate

        if 'results' in self.global_config:
            self.results = self.global_config['results']
        else:
            self.results = False

        ########################################################################
        ##                           QTL1 setting -- 1                        ##
        ########################################################################

        self.qtl1_type = self.global_config['input']['qtl1']['qtl_type']
        self.qtl1_file = self.global_config['input']['qtl1']['file']
        self.qtl1_biological_context = self.global_config['input']['qtl1']['biological_context'] 
        self._qtl1_type = self.global_config['input']['qtl1'].get('type', 'quant')

        if ('qtl_preprocessed_dir' in self.global_config['input']['qtl1']) and \
            self.global_config['input']['qtl1']['qtl_preprocessed_dir']:
            # 1. use custom eqtl preprocessed file
            qtl1_custom_preprocess_qtl_path = self.global_config['input']['qtl1']['qtl_preprocessed_dir']
        else:
            # 2. set eqtl output path
            qtl1_custom_preprocess_qtl_path = self.output_preprocessed_dir
            # e.g. preprocessed/eqtl/Whole_Blood/grouped
        self.qtl1_grouped_dir = os.path.join(qtl1_custom_preprocess_qtl_path, 
                                            self.qtl1_type, 
                                            self.qtl1_biological_context, 
                                            'grouped') 
        
        # e.g. preprocessed/eqtl/Whole_Blood/filtered_gene.tsv.gz
        self.qtl1_output_report = os.path.join(qtl1_custom_preprocess_qtl_path, 
                                            self.qtl1_type, 
                                            self.qtl1_biological_context, 
                                            'filtered_gene.tsv.gz')
        self.qtl1_LD_window = os.path.join(qtl1_custom_preprocess_qtl_path, 
                                            self.qtl1_type, 
                                            self.qtl1_biological_context, 
                                            'LD_window_per_gene.tsv.gz')
        self.qtl1_preprocesed_dir = os.path.join(qtl1_custom_preprocess_qtl_path,
                                                self.qtl1_type, 
                                            self.qtl1_biological_context)
        self.qtl1_col_dict = self.global_config['input']['qtl1']['col_name_mapping']

        self.qtl1_p_threshold = self.global_config['input']['qtl1']['p-value_threshold']
        if self.qtl1_p_threshold <= 0:
            self.qtl1_p_threshold = 1.0E-5



        self.qtl1_sep = self.global_config['input']['qtl1'].get('sep', '\t')
        if self.qtl1_sep == '\\t':
            self.qtl1_sep = '\t'
        elif (len(self.qtl1_sep)) > 1:
            logging.warning(f'QTL1 file separator is too long (can only be one char), using tab instead')
            self.qtl1_sep = '\t'



        self.qtl1_sample_size = self.global_config['input']['qtl1']['sample_size']


        ####################################################################
        ##                          QTL1 setting -- 2                     ##
        ####################################################################

        self.qtl2_type = self.global_config['input']['qtl2']['qtl_type']
        self.qtl2_file = self.global_config['input']['qtl2']['file']
        self.qtl2_biological_context = self.global_config['input']['qtl2']['biological_context'] 
        self._qtl2_type = self.global_config['input']['qtl2'].get('type', 'quant')

        if ('qtl_preprocessed_dir' in self.global_config['input']['qtl2']) and \
            self.global_config['input']['qtl2']['qtl_preprocessed_dir']:
            # 1. use custom eqtl preprocessed file
            qtl2_custom_preprocess_qtl_path = self.global_config['input']['qtl2']['qtl_preprocessed_dir']
        else:
            # 2. set eqtl output path
            qtl2_custom_preprocess_qtl_path = self.output_preprocessed_dir
            # e.g. preprocessed/eqtl/Whole_Blood/grouped
        self.qtl2_grouped_dir = os.path.join(qtl2_custom_preprocess_qtl_path, 
                                            self.qtl2_type, 
                                            self.qtl2_biological_context, 
                                            'grouped') 
        
        # e.g. preprocessed/eqtl/Whole_Blood/filtered_gene.tsv.gz
        self.qtl2_output_report = os.path.join(qtl2_custom_preprocess_qtl_path, 
                                            self.qtl2_type, 
                                            self.qtl2_biological_context, 
                                            'filtered_gene.tsv.gz')
        self.qtl2_LD_window = os.path.join(qtl2_custom_preprocess_qtl_path, 
                                            self.qtl2_type, 
                                            self.qtl2_biological_context, 
                                            'LD_window_per_gene.tsv.gz')
        self.qtl2_preprocesed_dir = os.path.join(qtl2_custom_preprocess_qtl_path,
                                                self.qtl2_type, 
                                            self.qtl2_biological_context)
        self.qtl2_col_dict = self.global_config['input']['qtl2']['col_name_mapping']

        self.qtl2_p_threshold = self.global_config['input']['qtl2']['p-value_threshold']
        if self.qtl2_p_threshold <= 0:
            self.qtl2_p_threshold = 1.0E-5



        self.qtl2_sep = self.global_config['input']['qtl2'].get('sep', '\t')
        if self.qtl2_sep == '\\t':
            self.qtl2_sep = '\t'
        elif (len(self.qtl2_sep)) > 1:
            logging.warning(f'QTL2 file separator is too long (can only be one char), using tab instead')
            self.qtl2_sep = '\t'

        self.qtl2_sample_size = self.global_config['input']['qtl2']['sample_size']



        self.tool_parent_dir = os.path.join(self.study_dir, 
                                            str(self.gwas_trait), 
                                            f"{self.qtl1_type}_{self.qtl1_biological_context}_{self.qtl2_type}_{self.qtl2_biological_context}_{self.population}",
                                            
                                            )
        
        self.coloc_base_dir = os.path.join(self.tool_parent_dir, 'coloc')
        self.susie_base_dir = os.path.join(self.tool_parent_dir, 'susie')
        self.rank_dir = os.path.join(self.tool_parent_dir, 'rank')
