source ~/.bashrc

raw_dir=../out/gwas_maxabsz/raw
processed_dir=../out/gwas_maxabsz/processed
mkdir -p ${raw_dir} ${processed_dir}

if [ ! -f ${raw_dir}/NCBI37.3.gene.loc ]; then
    wget -O ${raw_dir}/NCBI37.3.zip https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip 
    cd ${raw_dir}
    unzip NCBI37.3.zip && rm NCBI37.3.zip README REPORT
fi

python utils.py process_gwas_maxabsz_cli --raw_dir ${raw_dir} --processed_dir ${processed_dir}
