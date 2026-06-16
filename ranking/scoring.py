import os
from pathlib import Path

import pandas as pd
import scipy.stats as ss
import sys
import ranking.intact as intact
import ranking.rra as rra
import ranking.integrate as integrate
import numpy as np

def run_ranking(rpt_obj=None, output_file_path=None, out_integrate_results = None, prior_fun=None, sample_size=None, qtl_type=None):
    


    print(f"scoring {qtl_type}")
    out_dir = os.path.dirname(output_file_path) #/home/users/nus/e1124850/scratch/lc3test/processed/CONF/lc3test/Whole_Blood/EUR/sqtl/rank
    if not os.path.exists(out_dir):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
    geo_output_file = os.path.join(out_dir, 'geo.tsv') 
    intact_output_file = os.path.join(out_dir, 'intact.tsv')
    integrate_result = integrate.run(output_file_path=out_integrate_results, rpts=rpt_obj, qtl_type=qtl_type)
    geo_result = rra.run_ranking(output_file_path=geo_output_file, rpt=rpt_obj, sample_size=sample_size, method='GEO', qtl_type=qtl_type)
    intact_result = intact.run_ranking(rpt=rpt_obj, output_file_path=intact_output_file, prior_fun=prior_fun)

    geo_result_exist = geo_result is not None and os.path.exists(geo_result) and os.path.getsize(geo_result) > 0
    intact_result_exist = intact_result is not None and os.path.exists(intact_result) and os.path.getsize(
        intact_result) > 0

    if geo_result_exist and intact_result_exist:
        # output cols: ['phenotype_id', 'geo_ranking', 'avg_ranking', 'intact_probability']
        rgeo_df = pd.read_table(geo_result, usecols=['phenotype_id', 'geo_ranking'])
        intact_df = pd.read_table(intact_result,
                                usecols=lambda col: col in ['phenotype_id', 'avg_ranking', 'intact_probability'])
        result_df = pd.merge(left=rgeo_df, right=intact_df,
                            left_on='phenotype_id', right_on='phenotype_id',
                            how='outer')
        result_df['geo_ranking'] = ss.rankdata(result_df['geo_ranking'])
        result_df['geo_ranking'] = result_df['geo_ranking'].apply(lambda x: int(x))
        result_df['avg_ranking'] = ss.rankdata(result_df['avg_ranking'])
        result_df['avg_ranking'] = result_df['avg_ranking'].apply(lambda x: int(x))
    elif geo_result_exist:
        # intact_probability column is not present if user run TWAS only or colocalization methods
        result_df = pd.read_table(geo_result, usecols=['phenotype_id', 'geo_ranking'])
        result_df['geo_ranking'] = ss.rankdata(result_df['geo_ranking'])
        result_df['geo_ranking'] = result_df['geo_ranking'].apply(lambda x: int(x))

        result_df['avg_ranking'] = result_df['geo_ranking']
        os.remove(geo_result)
    elif intact_result_exist:
        result_df = pd.read_table(intact_result,
                                usecols=lambda col: col in ['phenotype_id', 'avg_ranking', 'intact_probability'])
        result_df['avg_ranking'] = ss.rankdata(result_df['avg_ranking'])
        result_df['avg_ranking'] = result_df['avg_ranking'].apply(lambda x: int(x))
        os.remove(intact_result)
    else:
        result_df = None


    if result_df is not None:
        result_df = result_df.round(4)
        if 'intact_probability' not in result_df.columns:
            result_df['intact_probability'] = np.nan

        result_df.to_csv(output_file_path, sep='\t', header=True, index=False, na_rep='NA')                
    return output_file_path