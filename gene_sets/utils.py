from urllib.request import urlopen
import pandas as pd
from bs4 import BeautifulSoup as bs
from os.path import join
import time
import fire
from tqdm import tqdm
import os
import numpy as np

from statsmodels.stats.multitest import multipletests
import scipy.stats

# utility function
def parse_twas_hub(prefix, url_base="http://twas-hub.org/traits"):
    
    # TODO: extract header information
    url = f"{url_base}/{prefix}"
    soup = bs(urlopen(url).read().decode("utf-8"), features='html.parser')
    trait_id = soup.h1.attrs['id']
    trait_name = soup.h1.contents[0]
    
    data = []
    table = soup.find('table', attrs={'id':'loci'})
    
    if table is None:
        return None
    
    # table head
    head = table.thead
    head = [ele.text.strip() for ele in head.find_all('th')]

    # table body
    body = table.find('tbody')
    rows = body.find_all('tr')
    for row in rows:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        data.append([ele for ele in cols if ele]) # Get rid of empty values
        
    df = pd.DataFrame(data, columns=head)
    return {'trait_id': trait_id, 'trait_name': trait_name, 'df': df}


# code for running depict
def run_depict(analysis_label, path_gwas_locus, results_dir):
    # START OF USER-SPECIFIC PARAMETERS 
    param_analysislabel = analysis_label
    
    # extract snp information for gwas locus
    with open(path_gwas_locus) as f:
        snps = f.readlines()[1].split(':')[1].strip().split(',')
        snps = [snp.strip() for snp in snps]

    path_snpfile = f'{results_dir}/{analysis_label}.snp.txt'
    with open(path_snpfile, 'w') as f:
        f.writelines('\n'.join(snps))

    path_results = results_dir
    flag_loci = 1
    flag_genes = 1
    flag_genesets = 0
    flag_tissues = 1
    param_ncores = 2
    # END OF USER-SPECIFIC PARAMETERS 

    path_locusgenerator_jar = "LocusGenerator/LocusGenerator.jar"
    path_depict_jar = "Depict/Depict.jar"

    # Construct loci
    if flag_loci:
        os.system(f"java -Xms512m -Xmx4000m -jar {path_locusgenerator_jar} LocusGenerator/config.xml {path_snpfile} {path_results}/{param_analysislabel} &> /dev/null")

    # Gene prioritization and reconstituted gene set enrichment analysis
    if flag_genes or flag_genesets:
        os.system(f"java -Xms512m -Xmx16000m -jar {path_depict_jar} {param_analysislabel} {flag_genes} {flag_genesets} 0 {param_ncores} {path_results} >> {path_results}/{param_analysislabel}.log 2>&1") 

    # Tissue/cell type enrichment analysis
    if flag_tissues:
        os.system(f"java -Xms512m -Xmx16000m -jar {path_depict_jar} {param_analysislabel} 0 1 1 {param_ncores} {path_results} >> {path_results}/{param_analysislabel}.log 2>&1") 

        
# processing code for different gene sets

# def process_twas(raw_dir, processed_dir):
#     """
#     Script for process TWAS gene sets. 
#     Download and parse all the TWAS gene sets from http://twas-hub.org/
#     """
#     # get all the prefix
#     url_base = "http://twas-hub.org/traits"
#     soup = bs(urlopen(url_base).read().decode("utf-8"), features='html.parser')
#     rows = soup.table.tbody.find_all('tr')
#     prefix_list = [row.find_all('td')[1].a.attrs['href'].split('/')[-1] for row in rows]
    
#     meta_info = {'PREFIX': [], 'ID': [], 'NAME': []}
#     for prefix in tqdm(prefix_list):
#         dic_twas = parse_twas_hub(prefix)
#         if dic_twas is None:
#             continue
#         dic_twas['df'].to_csv(join(raw_dir, f'{prefix}.csv'), index=False)
#         meta_info['PREFIX'].append(prefix)
#         meta_info['ID'].append(dic_twas['trait_id'])
#         meta_info['NAME'].append(dic_twas['trait_name'])
#         time.sleep(2)
        
#         # process
#         gene_sets = []
#         for i, row in dic_twas['df'].iterrows():
#             joint_genes = row['joint genes'].split(' ')
#             for gene in joint_genes:
#                 gene_sets.append([gene, row['best TWAS P']])
#         gene_sets = pd.DataFrame(gene_sets, columns=['GENE', 'SCORE'])
#         gene_sets.to_csv(join(processed_dir, f'{prefix}.csv'), index=False)
        
#     meta_info = pd.DataFrame(meta_info)
#     meta_info.to_csv(join(raw_dir, 'meta_info.csv'), index=False)

def process_twas(raw_dir, processed_dir):
    """
    Script for process TWAS gene sets. Raw data compiled by Kushal
    """
    with open(join(raw_dir, 'MAGMA_TWAS_matched_traits.txt')) as f:
        names = [line.strip().split() for line in f.readlines()]
        # convert from TWAS names to standard names (MAGMA)
        dic_name = {n[1] : n[0] for n in names}
    df = pd.read_csv(join(raw_dir, 'Genes_by_X_TWAS_abs_Z.txt'), sep='\t', index_col=0)
    df = df[list(dic_name.keys())].rename(columns=dic_name)
    
    for trait in df.columns:
        trait_df = df[[trait]].copy()
        trait_df.columns = ['z_score']
        trait_df['p_val'] = 2 * (1 - scipy.stats.norm.cdf(trait_df['z_score'].values))
        trait_df['fdr'] = multipletests(trait_df['p_val'].values, method='fdr_bh')[1]
        trait_df = pd.DataFrame({'GENE': trait_df.index, 'FDR': trait_df.fdr})
        trait_df.to_csv(join(processed_dir, f'{trait}.csv'), index=False)
    
    
def process_silver(raw_dir, processed_dir, subset):
    """
    Script for process causal gene sets as described in https://www.biorxiv.org/content/10.1101/2020.06.28.171561v1
    """
    assert subset in ['omim', 'drug']
    if subset == 'omim':
        df = pd.read_excel(join(raw_dir, 'table.xlsx'), sheet_name='S3', skiprows=2)
    elif subset == 'drug':
        df = pd.read_excel(join(raw_dir, 'table.xlsx'), sheet_name='S4', skiprows=3)
    
    
    df = df[['Trait', 'Gene Symbol']]
    traits = df['Trait'].unique()

    for trait in traits:
        trait_df = df[df['Trait'] == trait][['Gene Symbol']].drop_duplicates().reset_index(drop=True)
        trait_df.columns = ['GENE']
        trait_df['SCORE'] = 1
        trait_df.to_csv(join(processed_dir, f'{trait}.csv'), index=False)

def process_magma(raw_dir, processed_dir):
    
    df = pd.read_csv(join(raw_dir, 'MAGMA_v108_GENE_10_ZSTAT.txt'), sep='\t')
    for trait in df.columns:
        trait_df = df[[trait]].copy()
        trait_df.columns = ['z_score']
        # TODO: discuss with Martin if it make sense
        trait_df = trait_df[trait_df['z_score'] != 0]
        trait_df['p_val'] = 1 - scipy.stats.norm.cdf(trait_df['z_score'].values)
        trait_df['fdr'] = multipletests(trait_df['p_val'].values, method='fdr_bh')[1]
        trait_df = pd.DataFrame({'GENE': trait_df.index, 'FDR': trait_df.fdr})
        trait_df.to_csv(join(processed_dir, f'{trait}.csv'), index=False)

        
def process_depict(raw_dir, processed_dir):
    raw_dir = join(raw_dir, 'results_gwascatalog/')
    trait_list = [f[:f.rfind('_')] for f in os.listdir(raw_dir) if f.endswith('_geneprioritization.txt')]

    for trait in trait_list:
        df = pd.read_csv(join(raw_dir, f'{trait}_geneprioritization.txt'), sep='\t')
        df = df[['Gene symbol', 'Nominal P value']].rename(columns={'Gene symbol': "GENE", 'Nominal P value': "SCORE"})
        df['SCORE'] = multipletests(df['SCORE'].values, method='fdr_bh')[1]
        df.to_csv(join(processed_dir, f'{trait}.csv'), index=False)

        
def process_pops(raw_dir, processed_dir, subset):
    """
        subset: twas_hits / pops_hits / pops_all
    """
    assert subset in ['twas_hits', 'pops_hits', 'pops_all']
    if '_hits' in subset:
        ukb_df = pd.read_csv(join(raw_dir, 'UKB_AllMethods_GenePrioritization.txt.gz'), sep='\t')
        pass_df = pd.read_csv(join(raw_dir, 'PASS_AllMethods_GenePrioritization.txt.gz'), sep='\t')
        if subset == 'twas_hits':
            subset_score_col = 'twas_p'
            subset_rank_col = 'twas_rank'
        elif subset == 'pops_hits':
            subset_score_col = 'pops_score'
            subset_rank_col = 'pops_rank'
        for name, df in zip(['UKB', 'PASS'], [ukb_df, pass_df]):
            for trait in df['trait'].unique():
                trait_df = df[df['trait'] == trait].reset_index(drop=True)
                trait_df = trait_df[['gene', subset_score_col, subset_rank_col]].rename(
                    columns={'gene': 'GENE', 
                             subset_score_col: subset_score_col.upper(), 
                             subset_rank_col: subset_rank_col.upper()})
                trait_df.to_csv(join(processed_dir, f'{name}_{trait}.csv'), index=False)
                
    elif subset == 'pops_all':
        df = pd.read_csv(join(raw_dir, 'PoPS_FullResults.txt.gz'), sep='\t')
        df['cohort_trait'] = df['cohort'] + '_' + df['trait']
        for cohort_trait in df['cohort_trait'].unique():
            trait_df = df[df['cohort_trait'] == cohort_trait].reset_index(drop=True)
            trait_df = trait_df[['gene', 'pops_score']].rename(columns={'gene': 'GENE', 'pops_score': 'POPS_SCORE'})
            trait_df.to_csv(join(processed_dir, f'{cohort_trait}.csv'), index=False)
    else:
        raise NotImplementedError


def process_gene_max_abs_z(raw_dir, processed_dir, trait, window_size=10_000):
    sumstats_dir = '/n/holystore01/LABS/price_lab/Lab/ldsc/sumstats_formatted'
    reference_template = '/n/holystore01/LABS/price_lab/Lab/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.{}.bim'
    gene_loc_path = '/n/home12/khou/scTRS/gene_sets/out/raw/gwas_derived/NCBI37.3.gene.loc'
    
    # read reference
    reference_list = []
    for chr_i in range(1, 23):
        reference_list.append(pd.read_csv(reference_template.format(chr_i), sep='\t', usecols=[0, 1, 3], header=None, names=['CHR', 'SNP', 'BP']))
    reference = pd.concat(reference_list).reset_index(drop=True)
    
    # read gene location
    gene_loc = pd.read_csv(gene_loc_path, delim_whitespace=True, header=None, usecols=[1,2,3,5], names=['CHR', 'START', 'STOP', 'GENE'])
    gene_loc = gene_loc[gene_loc['CHR'].isin(np.arange(1, 23).astype(str))]
    gene_loc['CHR'] = gene_loc['CHR'].astype(int)
    
    def extract_gene_info(sumstats, gene, window_size):
        gene_index = (sumstats.CHR == gene.CHR) & (gene.START - window_size <= sumstats.BP) & (sumstats.BP < gene.STOP + window_size)
        gene_sumstats = sumstats[gene_index]
        return gene_sumstats['Z'].abs().max()
    
    sumstats = pd.read_csv(join(sumstats_dir, f'{trait}.sumstats'), sep='\t')
    sumstats = pd.merge(sumstats, reference, on='SNP')
    
    gene_loc['MAX_ABS_Z'] = gene_loc.apply(lambda gene : extract_gene_info(sumstats, gene, window_size), axis=1)
    
    gene_loc[['GENE', 'MAX_ABS_Z']].to_csv(join(processed_dir, f'{trait}.csv'), na_rep='NaN', index=False)


def process_hess(raw_dir, processed_dir):
    
    trait_list = [f[:-4] for f in os.listdir(raw_dir) if f.endswith('.csv')]

    for trait in trait_list:
        df = pd.read_csv(join(raw_dir, f'{trait}.csv'))
        df = df[['GENE', 'ESTIMATE']]
        df.to_csv(join(processed_dir, f'{trait}.csv'), index=False, na_rep='NaN')

    
if __name__ == '__main__':
    fire.Fire()