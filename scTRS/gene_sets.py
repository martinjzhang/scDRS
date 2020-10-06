import pandas as pd
import numpy as np
import os
from os.path import join
import glob
from gprofiler import GProfiler


class GeneSet():
    def __init__(self, root_dir, geneset_name, convert_mmusculus=True, verbose=False, **kwargs):
        """
        root_dir: Root directory to the gene set folder. See the preprocessing code for how to set this up.
        geneset_name: Name of the gene set
        trait: name of the traits
        kwargs: Parameters specific to different gene set, described as follows:
            Silver: 
        """
        assert geneset_name in ['depict', 
                                'silver_drug', 'silver_omim', 
                                'twas_kushal', 'magma_107', 'magma_108',
                                'pops_pops_hits', 'pops_twas_hits', 'pops_pops_all', 'gwas_max_abs_z', 'kegg', 'hess']
        
        # setting up the path
        self.data_dir = join(root_dir, 'processed', geneset_name)
        self.trait_list = [f[:-4] for f in os.listdir(self.data_dir) if f.endswith(".csv")]
        self.kwargs = kwargs
        self.verbose = verbose
        self.convert_mmusculus = convert_mmusculus
        self.geneset_name = geneset_name
        if verbose:
            print(f'Reading {len(self.trait_list)} from {geneset_name}')
            
        if self.convert_mmusculus:
            df_orth = pd.read_csv(join(root_dir, 'meta_info', 'hsapiens_mmusculus_gene_conversion.txt'), sep='\t')
            # Use only human genes that can be uniquely mapped to a mouse gene
            temp_df = df_orth.groupby(['incoming']).agg({'name':len})
            gene_list = list(temp_df.index[temp_df['name']==1])
            gene_list.sort()
            df_orth.index = df_orth['incoming']
            df_orth = df_orth.loc[gene_list]
            self.dic_convert_mmusculus = {df_orth['incoming'].values[x]:df_orth['name'].values[x] for x in np.arange(df_orth.shape[0])}

        
    def available_traits(self):
        """
            Get all the available traits
        """
        return self.trait_list
    
    def query_trait(self, trait):
        """
        Query the gene list for specific trait
        """
        
            
        gene_df = pd.read_csv(join(self.data_dir, f'{trait}.csv'))
        raw_len = gene_df.shape[0]

        if self.convert_mmusculus:
            gene_df['GENE'] = gene_df['GENE'].map(self.dic_convert_mmusculus)
            gene_df = gene_df.dropna()
        
        if self.verbose:
            print(f'#query_trait: {self.geneset_name}/{trait}, raw length: {raw_len}, after conversion to mmusculus: {gene_df.shape[0]}')
                    
        return gene_df


        