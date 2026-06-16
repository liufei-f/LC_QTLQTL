import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
import argparse

color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]

def color_map(LD):
    color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]
    if LD < 0.2:
        return color[0]
    elif LD < 0.4:
        return color[1]
    elif LD < 0.6:
        return color[2]
    elif LD < 0.8:
        return color[3]
    else:
        return color[4]

def matchsnp_preprocess(trait1_df, trait2_df, vcfdir, population, trait1_col_dict, trait2_col_dict):
    matchsnp = pd.concat([trait1_df[[trait1_col_dict['chrom'], 
                                     trait1_col_dict['snp'], 
                                     trait1_col_dict['pvalue'], 
                                     trait1_col_dict['beta'], 
                                     trait1_col_dict['position']]],
                         trait2_df[[trait2_col_dict['pvalue'],
                                   trait2_col_dict['beta']]]],
                         join='inner',axis=1)
    matchsnp.columns = ['chr', 'snp', 'trait1_p','trait1_beta','position','trait2_p','trait2_beta']
    matchsnp = matchsnp.sort_values('trait1_p')
    lead_snp = matchsnp.iloc[0,:]['snp']

    return lead_snp


def calculate_LD(trait_df, lead_snp, vcfdir, population, trait_col_dict):

    trait_df = trait_df[[trait_col_dict['chrom'],trait_col_dict['snp'],trait_col_dict['pvalue'], trait_col_dict['beta'], trait_col_dict['position']]]
    chrom = trait_df.iloc[0,:][trait_col_dict['chrom']]
    shell_command_plink_execute = 'plink --silent --vcf {} --r2 --ld-snp {} --ld-window 99999 --ld-window-kb 500 --ld-window-r2 0 --out {}'
    vcf_output_dir = f'{vcfdir}/{population}/chr{chrom}.vcf.gz'
    output_ld_file = f'/home/users/nus/e1124850/scratch/tmp/test'
    os.system(shell_command_plink_execute.format(vcf_output_dir, lead_snp, output_ld_file))
    ld = pd.read_csv(f'{output_ld_file}.ld',sep='\s+')
    ld.index = ld['SNP_B']
    ld = ld.drop_duplicates('SNP_B')
    trait_df["R2"] = trait_df.index.map(ld["R2"]).fillna(0)
    trait_df.columns = ['chr', 'snp', 'pvalue','beta','position','LD']
    trait_df = trait_df.fillna(0)
    trait_df['LD'][lead_snp] = 1
    trait_df['LD_color'] = trait_df['LD'].apply(lambda x: color_map(x))

    trait_df = trait_df.sort_values('position')
    trait_df['position'] = trait_df['position']/1e6

    return trait_df



def locuszoomplot(gene, chrom, vcfdir, population, 
                  trait1_type, trait1_name, trait1_path, trait1_col_dict, trait1_df,
                  trait2_type, trait2_name, trait2_dir, trait2_col_dict, 
                  outdir):
    # try:
        # trait1_df = pd.read_csv(f'{trait1_path}', sep='\t')
        # trait1_df.rename({v: k for k, v in trait2_col_dict.items()}, axis='columns', inplace=True)
        trait1_df.index = trait1_df[trait1_col_dict['snp']]
        
        trait1_df = trait1_df.sort_values(trait1_col_dict['pvalue'])
        if 'phenotype_id' in trait1_col_dict and 'phenotype_name' in trait1_col_dict:
            trait1_df = trait1_df[trait1_df[trait1_col_dict['phenotype_id']] == trait1_col_dict['phenotype_name']]
        
        trait1_df = trait1_df.drop_duplicates(trait1_col_dict['snp'])
        
    
        color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]
        classes = ['<0.2','0.2-0.4','0.4-0.6','0.6-0.8','>0.8']
        colors = ListedColormap(color)
        gene_ = gene.split('.')[0]

        trait2_df = pd.read_csv(\
            f'{trait2_dir}/{chrom}/{gene}.tsv.gz',
            sep='\t')
        # trait2_df.rename({v: k for k, v in trait2_col_dict.items()}, axis='columns', inplace=True)
        trait2_df = trait2_df.sort_values(trait2_col_dict['pvalue'])
        trait2_df.index = trait2_df[trait2_col_dict['snp']]
        trait2_df = trait2_df.drop_duplicates(trait2_col_dict['snp'])

        # ld_window_ls = ld_window.split(',')
        # fixed_window_ls = fixed_window.split(',')
        lead_snp = matchsnp_preprocess(trait1_df, trait2_df, vcfdir, population, trait1_col_dict, trait2_col_dict)


        trait1_df = calculate_LD(trait1_df, lead_snp, vcfdir, population, trait1_col_dict)
        trait2_df = calculate_LD(trait2_df, lead_snp, vcfdir, population, trait2_col_dict)

        pos_min = min(min(trait1_df['position']), min(trait2_df['position']))
        pos_max = max(max(trait1_df['position']), max(trait2_df['position']))


        trait1_df['pvalue'] = -np.log10(trait1_df['pvalue'])
        trait2_df['pvalue'] = -np.log10(trait2_df['pvalue'])
        ################################################
        #                  1. locuszoom                #
        ################################################
        # trait2
        fig = plt.figure(figsize = (5.5, 4.5))
        ax1 = plt.subplot(211)
        ax1.scatter(trait1_df['position'],trait1_df['pvalue'], 
                    c=trait1_df['LD_color'], cmap=colors, s=50, 
                    edgecolors='grey', linewidths=0.01)
        color_ = 'darkviolet'
        ax1.scatter(trait1_df.loc[lead_snp]['position'], 
                trait1_df.loc[lead_snp]['pvalue'], c='darkviolet', 
                marker='D', s=60, edgecolors='black', linewidths=1)
        ax1.text(trait1_df.loc[lead_snp]['position'], 
                trait1_df.loc[lead_snp]['pvalue'], f"{lead_snp}", fontsize=10)
        ax1.set_xlim(pos_min, pos_max)
        # ## trait2 ld-based loci test region
        # ax1.axvline(x=ld_window_ls[0]/1e6, color='cornflowerblue', linestyle='--', linewidth=1.5)
        # ax1.axvline(x=ld_window_ls[1]/1e6, color='cornflowerblue', linestyle='--', linewidth=1.5)
        # ## global ld-based loci test region
        # ax1.axvline(x=fixed_window_ls[0]/1e6, color='tan', linestyle='--', linewidth=1.5)
        # ax1.axvline(x=fixed_window_ls[1]/1e6, color='tan', linestyle='--', linewidth=1.5)
        

        ax1.set_ylabel(f'{trait1_type}\n{trait1_name} ' + '-log$_{10}$(${P}$)',fontsize=13)
        ax1.spines.right.set_visible(False)
        ax1.spines.top.set_visible(False)
    
        # locus trait1
        ax2 = plt.subplot(212)
        ax2.scatter(trait2_df['position'],trait2_df['pvalue'], 
                    c=trait2_df['LD_color'], cmap=colors, s=50, 
                    edgecolors='grey', linewidths=0.01)
        ax2.scatter(trait2_df.loc[lead_snp]['position'], 
                    trait2_df.loc[lead_snp]['pvalue'], c='darkviolet', 
                    marker='D',s=60, edgecolors='black', linewidths=1)
        ax2.text(trait2_df.loc[lead_snp]['position'], 
                trait2_df.loc[lead_snp]['pvalue'], f"{lead_snp}", fontsize=10)
        ax2.set_xlabel(f'Chromosome {chrom} (MB)',fontsize=13)
        ax2.set_xlim(pos_min, pos_max)
        # ## trait2 ld-based loci test region
        # ax2.axvline(x=ld_window_ls[0]/1e6, color='cornflowerblue', linestyle='--', linewidth=1.5)
        # ax2.axvline(x=ld_window_ls[1]/1e6, color='cornflowerblue', linestyle='--', linewidth=1.5)
        # ## global ld-based loci test region
        # ax2.axvline(x=fixed_window_ls[0]/1e6, color='tan', linestyle='--', linewidth=1.5)
        # ax2.axvline(x=fixed_window_ls[1]/1e6, color='tan', linestyle='--', linewidth=1.5)
        ax2.set_ylabel(f'{trait2_type}\n{trait2_name}'+'-log$_{10}$(${P}$)',fontsize=13)

        ax2.spines.right.set_visible(False)
        ax2.spines.top.set_visible(False)
        # plt.savefig(f'{outdir}/locuszoom_{trait1_type}_{trait1_col_dict['phenotype_name']}_{trait2_type}_{trait2_name}_{gene}.pdf',
        #             format='pdf',bbox_inches='tight')
        plt.savefig(f'{outdir}/locuszoom_{trait1_type}_{trait1_col_dict['phenotype_name']}_{trait2_type}_{trait2_name}_{gene}.png',
                    format='png', dpi=300, bbox_inches='tight')
        plt.close()
        
    # except Exception as e:
    #     print(f"Error: {e}")
    #     print(f"Gene: {gene}, chrom: {chrom}")


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--gene', dest='gene',
#                         help='')
#     parser.add_argument('--chrom', dest='chrom',
#                         help='')
#     parser.add_argument('--vcfdir', dest='vcfdir',
#                         help='')
#     parser.add_argument('--population', dest='population',
#                         help='')

#     parser.add_argument('--trait2_trait', dest='trait2_trait',
#                        help='')
#     parser.add_argument('--trait2_path', dest='trait2_trait',
#                        help='')
#     parser.add_argument('--trait1_dir', dest='trait2_trait',
#                        help='')
#     parser.add_argument('--ld_window', dest='ld_window', default='59428243,59530853',
#                        help='')
#     parser.add_argument('--fixed_window', dest='fixed_window',default='59428243,59530853',
#                        help='')
#     parser.add_argument('--outdir', dest='outdir',
#                        help='')


#     args = parser.parse_args()
#     print(f'Accepted args:\n {args}')

#     locuszoomplot(args.gene, args.chrom, args.vcfdir, args.population, args.trait2_trait, args.trait2_path, args.trait1_dir, args.ld_window, args.fixed_window, args.outdir)



vcfdir = '/home/users/nus/e1124850/e1124850/locuscompare3/1000genomes_newrsid'
population = 'EUR'

trait1_type = 'caQTL'
trait1_name = 'multitissue'
trait1_path = '/home/users/nus/e1124850/scratch/caqtl/global_caQTL_mapping_data/preprocessed_wenz_2025_caqtl.tsv.gz'

trait1_df = pd.read_csv(f'{trait1_path}', sep='\t')

trait1_col_dict = {'chrom': 'chr', 
                 'position': 'pos', 
                 'alt': 'alt', 
                 'ref': 'ref', 
                 'snp': 'rsid', 
                 'beta':'slope',
                 'se':'slope_se',
                 'pvalue': 'pval_nominal',
                 'phenotype_id': 'peak_region',
                 'phenotype_name': 'chr22_22809848_22810933'}



chrom = '22'
trait2_type = 'eQTL'
trait2_name = 'Whole_Blood'
trait2_dir = f'/home/users/nus/e1124850/scratch/lc3testgwas/preprocessed/{trait2_type.lower()}/{trait2_name}/grouped/'
trait2_col_dict = {'chrom': 'chr', 
                 'position': 'pos', 
                 'alt': 'alt', 
                 'ref': 'ref', 
                 'snp': 'rs_id_dbSNP155_GRCh38p13', 
                 'beta':'slope',
                 'se':'slope_se',
                 'pvalue': 'pval_nominal'}


outdir = '/home/users/nus/e1124850/lc3private'

gene = 'ENSG00000211669.3'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

gene = 'ENSG00000278196.3'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

gene = 'ENSG00000211670.2'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

gene = 'ENSG00000211666.2'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

gene = 'ENSG00000211662.2'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

gene = 'ENSG00000211667.3'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


gene = 'ENSG00000211668.2'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

gene = 'ENSG00000288861.1'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

# df[df['phenotype_id_qtl2'] == 'chr22_22809848_22810933'][['chrom','phenotype_id','overall_H4']]
#          chrom                               phenotype_id  overall_H4
# 0           22  ENSG00000211669.3_chr22_22809848_22810933      1.0000
# 19          22  ENSG00000278196.3_chr22_22809848_22810933      1.0000
# 21          22  ENSG00000211670.2_chr22_22809848_22810933      1.0000
# 307431      22  ENSG00000211666.2_chr22_22809848_22810933      0.9322
# 883587      22  ENSG00000211662.2_chr22_22809848_22810933      0.0803
# 1167998     22  ENSG00000211667.3_chr22_22809848_22810933      0.0249
# 1523535     22  ENSG00000211668.2_chr22_22809848_22810933      0.0005
# 1606660     22  ENSG00000288861.1_chr22_22809848_22810933      0.0000


outdir = '/home/users/nus/e1124850/lc3private'

trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr20_32739914_32739979'}

gene = 'ENSG00000149600.12'
chrom = '20'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr20_32739220_32739604'}

gene = 'ENSG00000149600.12'
chrom = '20'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr8_6970364_6971852'}
gene = 'ENSG00000215378.3'
chrom = '8'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr8_7020538_7021366'}
gene = 'ENSG00000215378.3'
chrom = '8'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr8_7018081_7018184'}
gene = 'ENSG00000215378.3'
chrom = '8'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)



trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr8_6969562_6969658'}
gene = 'ENSG00000215378.3'
chrom = '8'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)



trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr8_6968852_6969243'}
gene = 'ENSG00000215378.3'
chrom = '8'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)



trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr7_75724488_75725296'}
gene = 'ENSG00000127946.17'
chrom = '7'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)



trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr1_161636282_161636678'}
gene = 'ENSG00000162747.12'
chrom = '1'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr1_161611179_161613398'}
gene = 'ENSG00000162747.12'
chrom = '1'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)



trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr3_128950108_128951251'}
gene = 'ENSG00000287110.2'
chrom = '3'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)




trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr4_184850280_184851727'}
gene = 'ENSG00000286256.2'
chrom = '4'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)



trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr10_826967_827367'}
gene = 'ENSG00000185736.16'
chrom = '10'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr8_7014872_7015090'}
gene = 'ENSG00000223629.1'
chrom = '8'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr11_4175617_4176340'}
gene = 'ENSG00000167325.15'
chrom = '11'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)


trait1_col_dict = {'chrom': 'chr', 'position': 'pos', 'alt': 'alt', 'ref': 'ref', 
                   'snp': 'rsid', 'beta':'slope','se':'slope_se','pvalue': 'pval_nominal',
                   'phenotype_id': 'peak_region',
                   'phenotype_name': 'chr3_128686462_128687368'}
gene = 'ENSG00000231305.4'
chrom = '3'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

gene = 'ENSG00000163902.12'
chrom = '3'
locuszoomplot(gene, chrom, vcfdir, population, trait1_type, trait1_name, 
              trait1_path, trait1_col_dict, trait1_df, trait2_type, trait2_name, trait2_dir, trait2_col_dict, outdir)

