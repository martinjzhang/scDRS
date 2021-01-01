source ~/.bashrc

raw_dir=../out/magma/raw
processed_dir=../out/magma/processed
mkdir -p ${raw_dir} ${processed_dir}

cp /n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/gene_annotation/MAGMA-v108/MAGMA_v108_GENE_10_PSTAT.txt ${raw_dir}
cp /n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/gene_annotation/MAGMA-v108/MAGMA_v108_GENE_10_ZSTAT.txt ${raw_dir}
cp /n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/sumstats/Description_080419.xlsx ${raw_dir}/description.xlsx

python utils.py process_magma --raw_dir ${raw_dir} --processed_dir ${processed_dir}
