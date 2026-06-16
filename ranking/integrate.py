import os
import pandas as pd
import logging
from pathlib import Path
from ranking.constants import PHENOTYPE_ID_COL_NAME, GENE_ID_COL_NAME
from ranking.constants import TOOL_SIG_COL_INFO
from ranking.constants import RESULT_TYPE_PVAL
import sys

RANK_COL_NAME = 'rank'

def run(output_file_path, rpts, qtl_type):
    print(f"Integrating results {qtl_type}")
    if rpts is None or len(rpts) == 0:
        logging.warning('No reports provided')
        return None
    if not any(rpts.values()):
        logging.warning('All reports are empty or none')
        return None
    df_list = []
    for tool, sig_column, sig_type, lead_variant in TOOL_SIG_COL_INFO:
        tool_rpt = rpts.get(tool)
        if tool_rpt is None or len(tool_rpt) == 0:
            continue
        df_list.append(read_tool_result(rpt = tool_rpt, tool_name=tool, sig_col_name=sig_column, sig_type=sig_type, lead_variant = lead_variant, qtl_type=qtl_type))

    integrate_df = None
    for rpt_df in [df for df in df_list if df is not None and df.shape[0] > 0]:
        if integrate_df is None:
            integrate_df = rpt_df
        else:
            integrate_df = pd.merge(left=integrate_df, right=rpt_df,
                                  left_on='lead_variant', right_on='lead_variant',
                                  how='outer')
    if integrate_df is None:
        return None

    integrate_df.to_csv(output_file_path, sep='\t', header=True, index=False, na_rep='NA')
    return output_file_path


def read_tool_result(rpt, tool_name, sig_col_name, sig_type, lead_variant, qtl_type):
    
    if rpt is None or (not os.path.exists(rpt)) or os.path.getsize(rpt) <= 0:
        return None
    # rpt_df = pd.read_table(rpt, usecols=[sig_col_name, PHENOTYPE_ID_COL_NAME])
    rpt_df = pd.read_table(rpt)
    print(f"[DEBUG] Columns in rpt_df: {rpt_df.columns.tolist()}")
    print(f"qtltype: {qtl_type}")

    if tool_name in ['coloc','ecaviar','fastenloc','smr']:
        rpt_df = rpt_df[[sig_col_name, PHENOTYPE_ID_COL_NAME,lead_variant]]
        # rpt_df['ix'] = rpt_df[PHENOTYPE_ID_COL_NAME] + '_' + rpt_df[lead_variant]
        rpt_df['ix'] = rpt_df[PHENOTYPE_ID_COL_NAME] + '_' + rpt_df[lead_variant]
        rpt_df.drop_duplicates(subset=['ix'], inplace=True)
        rpt_df.sort_values(sig_col_name, ascending=sig_type == RESULT_TYPE_PVAL, inplace=True)
        rpt_df[tool_name] = rpt_df[sig_col_name]
        rpt_df['qtl_type'] = qtl_type
        rpt_df.drop(labels=[sig_col_name], axis=1, inplace=True)
        return rpt_df

    else:
        rpt_df = rpt_df[[sig_col_name, GENE_ID_COL_NAME, PHENOTYPE_ID_COL_NAME]]
        rpt_df['ix'] = rpt_df[PHENOTYPE_ID_COL_NAME] ## ??
        rpt_df.drop_duplicates(subset=['ix'], inplace=True)
        rpt_df.sort_values(sig_col_name, ascending=sig_type == RESULT_TYPE_PVAL, inplace=True)
        rpt_df[tool_name] = rpt_df[sig_col_name]
        rpt_df.drop(labels=[sig_col_name], axis=1, inplace=True)
        return rpt_df