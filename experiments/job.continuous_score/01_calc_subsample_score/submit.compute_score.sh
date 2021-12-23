for f in geneset/*.gs.batch; do
    prefix="$(basename -- $f)"
    prefix=${prefix::-9}
    echo $prefix
    n_batch=$(ls -l geneset/${prefix}.gs.batch | wc -l)
    sbatch --array 0-$((n_batch-2)) compute_score.sh ${prefix} vs
    sbatch --array 0-$((n_batch-2)) compute_score.sh ${prefix} uniform
done
