# README

```bash
n_batch=37
for dataset in ayhan_2021 habib_2016 habib_2017 yao_2021 zeisel_2015 zeisel_2018 zhong_2020; do
    sbatch --array 0-${n_batch} score_${dataset}.sh
done
```