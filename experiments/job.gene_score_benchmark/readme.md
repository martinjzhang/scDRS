# Gene score benchmark

## Make gene set files
```python
import pandas as pd

meta_df = pd.DataFrame({"method": ["gwas_maxabsz", "magma", "hess"],
                        "col": ["MAX_ABS_Z", "FDR", "ESTIMATE"],
                        "ascending": [False, True, False]})

for i, row in meta_df.iterrows():
    df_gs = pd.DataFrame(columns=['TRAIT', 'GENESET'])
    gene_score_path = f"/n/home12/khou/scTRS/gene_score/out/{row.method}/processed"
    for file in os.listdir(gene_score_path):
        trait = file.replace('.csv','')
        temp_df = pd.read_csv(join(gene_score_path, file), sep=',')
        temp_df = temp_df.loc[~temp_df[row.col].isna()]
        temp_df = temp_df.sort_values(by=row.col, ascending=row.ascending)
        df_gs.loc[trait] = [trait,','.join(temp_df['GENE'][0:500])]
        
    df_gs.to_csv(f"out/gs_file/{row.method}.top500.gs", sep='\t', index=False)
    BATCH_SIZE = 4
    for i_batch in range(np.ceil(df_gs.shape[0] / BATCH_SIZE).astype(int)):
        df_gs.iloc[i_batch*BATCH_SIZE:(i_batch+1)*BATCH_SIZE].to_csv(f"out/gs_file/{row.method}.top500.batch{i_batch}.gs",
                    sep='\t', index=False)
```

## Submit the jobs to get the score

### Step 1: submit all the scores at once
```bash
for method in gwas_maxabsz hess magma; do
    n_batch=$(ls -l out/gs_file/${method}.*.batch*.gs | wc -l)
    sbatch --array=0-$((n_batch - 1)) compute_score.sh ${method}
done
```

### Step 2
Jobs sometimes does not finish. We submit the unfished ones with the following code.
```python
import os
import glob
import re
import subprocess
import pandas as pd

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

for method in ["gwas_maxabsz", "hess", "magma"]:
    gs_file_list = natural_sort(glob.glob(f"out/gs_file/{method}.*.batch*.gs"))
    n_batch = len(gs_file_list)
    todo_batch_list = []
    for batch_i in range(n_batch):
        gs_file = gs_file_list[batch_i]
        trait_list = pd.read_csv(gs_file, sep='\t')['TRAIT'].values
        for trait in trait_list:
            if not os.path.exists(f"out/score_file/score.facs.{method}.top500/{trait}.score.gz"):
                todo_batch_list.append(batch_i)
    todo_batch_list = sorted(list(set(todo_batch_list)))
    print(method, todo_batch_list)
    for batch_i in todo_batch_list:
        cmd = f"sbatch --array={batch_i} compute_score.sh {method}"
        print(cmd)
        subprocess.call(cmd, shell=True)   
```

## Analyze the results
See `gene_score_benchmark.ipynb` for the results analysis.