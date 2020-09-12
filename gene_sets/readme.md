# Disease associated gene sets
We compile several existing trait associated gene sets for integrating scRNA-seq and GWAS dataset.

Ideally, each gene set has information for 1. trait / disease 2. set of genes 3. strength of the association, e.g. p-value. Therefore, each gene set is a file with these three columns.

Below is a summary we obtain from different gene sets:
- TWAS: We download the summary statistics from http://twas-hub.org/, using the `best TWAS P` as the score.
- DEPICT: We ran the DEPICT software, and run on 63 GWAS Catalog downloaded from DEPICT website. The nominal p-values are converted to FDR as the scores.
- Silver gene sets: We use the gene sets compiled from https://www.biorxiv.org/content/10.1101/2020.06.28.171561v1. There are OMIM, which are known mendelian genes to the diseases, and known drug targets. Each gene are assigned equal score 1.
- MAGMA: Nominal p-values are converted to FDR, which is taken as the score.
- PoPS: the scores are directly used.
In the follows, we describe how we processed to get the gene sets.


## TWAS
Downloaded from http://twas-hub.org/. 
```bash
mkdir -p out/raw/twas out/processed/twas
python utils.py process_twas --raw_dir out/raw/twas --processed_dir out/processed/twas
```


## Silver genesets 
The gene sets is taken from https://www.biorxiv.org/content/10.1101/2020.06.28.171561v1
Taken from Table S3 (List of genes whose perturbation is known to cause a Mendelian form of the disease, or influence the trait) and Table S4 (List of drug targets for each disease or trait)
```
mkdir -p out/raw/silver
wget -O out/raw/silver/table.xlsx https://www.biorxiv.org/content/biorxiv/early/2020/06/28/2020.06.28.171561/DC2/embed/media-2.xlsx
for subset in drug omim; do
    mkdir -p out/processed/silver_${subset}
    python utils.py process_silver --raw_dir out/raw/silver --processed_dir out/processed/silver_${subset} --subset ${subset}
done

```

## DEPICT

```bash
mkdir -p out/depict/processed out/depict/raw
cd out/depict/raw
# download data
wget https://data.broadinstitute.org/mpg/depict/depict_140721.tar.bz2 
tar xvfj depict_140721.tar.bz2
wget https://data.broadinstitute.org/mpg/depict/gwascatalog_140201.tar.gz
tar xvf gwascatalog_140201.tar.gz
```
There are some files in `gwascatalog_140201` with parenteses, let's remove them.
```bash
cd gwascatalog_140201
for i in ./*
do
     name=`echo $i | sed 's/[,()]//g'`
     mv "$i" "$name"
done
```
Now copy `utils.py` to `out/depict/raw`

Use the following script, save this as `submit.sh` in `out/depict/raw`
```bash
trait=$1

results_dir=results_gwascatalog
mkdir -p ${results_dir}
python utils.py run_depict --analysis_label ${trait} --path_gwas_locus gwascatalog_140201/${trait}.txt --results_dir ${results_dir}
```

Now submit the jobs
```bash
file_list=$(ls gwascatalog_140201)
for f in ${file_list}; do
    trait="$(basename $f .txt)"
    echo $trait
    bash submit.sh ${trait}
done
```

Wait untils the jobs are done. 

```
python utils.py process_depict --raw_dir out/raw/depict --processed_dir out/processed/depict
```

## MAGMA
Starting with files provides by Kushal. We format the data for later convenience.
```bash
mkdir -p out/raw/magma out/processed/magma
cp /n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/gene_annotation/Genes_by_X_kushal/Genes_by_X_MAGMA_10kb_Z.txt out/raw/magma/magma_10kb_Z.txt
cp /n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/sumstats/Description_080419.xlsx out/raw/magma/description.xlsx

python utils.py process_magma --raw_dir out/raw/magma --processed_dir out/processed/magma
```

## PoPS

```bash
mkdir -p out/raw/pops out/processed/pops_pops out/processed/pops_twas
wget -O out/raw/pops/PoPS_FullResults.txt.gz https://www.dropbox.com/sh/dz4haeo48s34sex/AADyog193x9waw0YXgdxDobja/results/PoPS_FullResults.txt.gz
wget -O out/raw/pops/PASS_AllMethods_GenePrioritization.txt.gz https://www.dropbox.com/sh/dz4haeo48s34sex/AAARC3I2C3mMPl1wmLlFzY7ta/results/PASS_AllMethods_GenePrioritization.txt.gz
wget -O out/raw/pops/UKB_AllMethods_GenePrioritization.txt.gz https://www.dropbox.com/sh/dz4haeo48s34sex/AAAhP3P1f8RzIiTjgcArHEvIa/results/UKB_AllMethods_GenePrioritization.txt.gz

python utils.py process_pops --raw_dir out/raw/pops --processed_dir out/processed/pops_pops -subset pops
python utils.py process_pops --raw_dir out/raw/pops --processed_dir out/processed/pops_twas -subset twas

```

