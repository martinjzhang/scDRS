# Aggregating SNP-level association statistics into gene-level statistics.

We compiled a set of gene-level statistics for aggregating SNP-level association statistics into gene-level statistics and benchmark how different strategies impacts cell scoring procedure.

Below is a summary we obtain using different methods:
- MAGMA: v1.08, with 10kb window. Starting with files provides by Kushal. Run `bash magma.sh`
- GWAS maximum absolute Z-score: For each gene, we take a 10-kb window, take the maximum absolute z-score as the score for this gene. Run `bash gwas_maxabsz.sh`
- HESS: A gene-level heritability estimatation method.
