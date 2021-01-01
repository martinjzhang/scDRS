import pandas as pd
import numpy as np
import os
from os.path import join
import glob
from gprofiler import GProfiler


class GeneScore():
    def __init__(self, root_dir, method_name, convert_mmusculus=True, verbose=False, **kwargs):
        """
        root_dir: Root directory to the gene set folder. See the preprocessing code for how to set this up.
        method_name: Name of the gene set
        trait: name of the traits
        kwargs: Parameters specific to different gene set, described as follows:
            Silver: 
        """
        assert method_name in ['MAGMA_v108', 'MAGMA_v107', 'GWAS_MAX_ABS_Z', 'HESS']
        method_dir_lookup = {
            'MAGMA_v108': 'magma_108',
            'MAGMA_v107': 'magma_107',
            'GWAS_MAX_ABS_Z': 'gwas_max_abs_z',
            'HESS': 'hess'
        }
        
        # setting up the path
        self.data_dir = join(root_dir, 'processed', method_dir_lookup[method_name])
        self.trait_list = [f[:-4] for f in os.listdir(self.data_dir) if f.endswith(".csv")]
        self.kwargs = kwargs
        self.verbose = verbose
        self.convert_mmusculus = convert_mmusculus
        self.method_name = method_name
        if verbose:
            print(f'Reading {len(self.trait_list)} from {method_name}')
            
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
    
    def map_to_mmusculus(self, gene_list):
        """
        Given a list of genes is hsapiens, convert that to mmusculus
        """
        converted = [self.dic_convert_mmusculus[g] for g in gene_list if g in self.dic_convert_mmusculus]
        return converted
    
    def query_trait(self, trait, sort_by=None, ascending=None, num_genes=None):
        """
        Query the gene list for specific trait
        trait: trait name to query
        sorted: whether return a sorted data frame
        num_genes: number of top genes to retrieve
        """
        
            
        gene_df = pd.read_csv(join(self.data_dir, f'{trait}.csv'))
        raw_len = gene_df.shape[0]

        if self.convert_mmusculus:
            gene_df['GENE'] = gene_df['GENE'].map(self.dic_convert_mmusculus)
            gene_df = gene_df.dropna()
        
        if self.verbose:
            print(f'#query_trait: {self.method_name}/{trait}, raw length: {raw_len}, after conversion to mmusculus: {gene_df.shape[0]}')
        
        if sort_by is not None:
            gene_df = gene_df.sort_values(by=sort_by, ascending=ascending)
        
        if num_genes is not None:
            gene_df = gene_df.iloc[0 : num_genes, :].reset_index(drop=True)
            
        return gene_df


        