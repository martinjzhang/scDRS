Here is a guide to compute MAGMA gene sets directly from GWAS summary statistics:

```shell
# Step 1: download MAGMA software, gene location file, and reference data from
# https://ctg.cncr.nl/software/magma after this step, one should have a folder <MAGMA_DIR>
# with the following files:
# 1) <MAGMA_DIR>/magma 2) <MAGMA_DIR>/g1000_eur.(bed|bim|fam) 3) <MAGMA_DIR>/NCBI37.3.gene.loc

magma_dir=<MAGMA_DIR>

# Step 2: make gene annotation file for MAGMA using the following command, this only needs to be done
# once for different GWAS summary statistics, and the file will be saved to out/step1.genes.annot
mkdir -p out/step1
${magma_dir}/magma \
    --annotate window=10,10 \
    --snp-loc ${magma_dir}/g1000_eur.bim \
    --gene-loc ${magma_dir}/NCBI37.3.gene.loc \
    --out out/step1

# Step 3: run MAGMA using the following command, this takes a GWAS file ${trait}.pval,
# which at least has the columns: SNP, P, N, which corresponds to the SNP id
# (matched to the ${magma_dir}/g1000_eur.bim), p-value, sample size. For example,
# <trait>.pval file looks like
#
# CHR     BP      SNP             P           N
# 1       717587  rs144155419     0.453345    279949
# 1       740284  rs61770167      0.921906    282079
# 1       769223  rs60320384      0.059349    281744
#
# After this step, one should obtain a file out/step2/${trait}.gene.out, and the top genes with
# largest Z-scores can be input to scDRS.

mkdir -p out/step2
${magma_dir}/magma \
    --bfile ${magma_dir}/g1000_eur \
    --pval ${trait}.pval use='SNP,P' ncol='N' \
    --gene-annot out/step1.genes.annot \
    --out out/step2/${trait}
```

We can then map IDs in `out/step2/${trait}.gene.out` to gene symbols using the map from `<MAGMA_DIR>/NCBI37.3.gene.loc` 
