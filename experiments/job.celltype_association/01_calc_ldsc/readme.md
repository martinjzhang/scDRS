# README

Comparing LDSC-SEG with scTRS results for finding celltype - trait relationship.

## Run LDSC-SEG

### Prepare gene location
```bash
mkdir out/gene_loc
cd out/gene_loc
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip && unzip NCBI37.3.zip
```

```python
import pandas as pd
gene_loc = pd.read_csv("NCBI37.3.gene.loc", names=["ID", "CHR", "START", "END", "STRAND", "GENE"], delim_whitespace=True)
gene_loc = gene_loc[~gene_loc.CHR.isin(["X", "Y"])]
gene_loc["CHR"] = "chr" + gene_loc["CHR"]
gene_loc[["GENE", "CHR", "START", "END"]].to_csv("gene_loc.txt", sep='\t', index=False)
```

### Prepare the gene sets

### Calculate LD score for gene sets

```bash

for dset in ts_facs tms_facs tms_droplet ; do
    n_ct=$(cat out/${dset}/ct_list.txt | wc -l)
    sbatch --array=1-${n_ct} calc_ldscore.sh ${dset}
done
```

```bash
for dset in ts_facs tms_facs tms_droplet; do
    # to check if anything unfinishes
    n_ct=$(cat out/${dset}/ct_list.txt | wc -l)
    for i in $(seq 1 $n_ct); do
        ct=$(sed "${i}q;d" out/${dset}/ct_list.txt)
        if [ ! -f out/${dset}/ldscore/${ct}.22.l2.ldscore.gz ]; then
            sbatch --array=$i calc_ldscore.sh ${dset}
        fi
    done
done
```
### run LD score

```bash
n_trait=$(cat out/trait_list.txt | wc -l)
for dset in ts_facs tms_facs tms_droplet; do
    sbatch --array=1-${n_trait} calc_cts.sh ${dset}
done

for i in $(seq 1 $n_trait); do
    trait=$(sed "${i}q;d" out/trait_list.txt)
    if [ ! -f out/cts_result/${trait}.cell_type_results.txt ]; then
        sbatch -t $i calc_cts.sh
    fi
done
```

# Summarize results for full LDSC-SEG outputs
```python
import pandas as pd
from os.path import join
out_dir = "out/"

with open(join(out_dir, "trait_list.txt")) as f:
    trait_list = [trait.strip() for trait in f.readlines()]
    
with open(join(out_dir, "tms_facs.ldcts")) as f:
    ct_list = [line.split()[0] for line in f.readlines()]

for trait in trait_list:
    print(trait)
    df = []
    for ct in ct_list:
        rls = pd.read_csv(join(out_dir, "cts_result_full", trait, f"{ct}.results"), delim_whitespace=True).iloc[0]
        rls["Name"] = ct
        df.append(rls)
    df = pd.concat(df, axis=1).T.set_index("Name")
    df.to_csv(join(out_dir, "cts_result_full", f"{trait}.cell_type_results.txt"), sep='\t')
```
