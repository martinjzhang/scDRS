# Subset of code and data to reproduce main results of the paper

First download the data [scDRS_data_release_092121](XXX).  

### Compute scDRS scores for TMS FACS + 74 diseases/traits
- Score files were already included in `scDRS_data_release_092121/score_file/score.tms_facs_with_cov.magma_10kb_1000`.
- You can also compute them yourself by setting `DATA_PATH` to your local folder of `scDRS_data_release_092121` and run `job.reproduce/reproduce_compute_score.tms_facs_with_cov.magma_10kb_1000.sh`.

### Cell type-level analysis (Fig. 3)
XXX

### T cell analysis (Fig. 4): 
- Fig. 4A,B,C: set `DATA_PATH` to your local folder of `scDRS_data_release_092121` and run `job.reproduce/reproduce_tcell.ipynb`.
- Fig. 4D: set `DATA_PATH` to your local folder of `scDRS_data_release_092121` and run `job.reproduce/reproduce_tcell_gene.ipynb`. 

### Neuron analysis (Fig. 5A,B):
XXX

### Hepatocyte analysis (Fig. 5C,D): 
Set `DATA_PATH` to your local folder of `scDRS_data_release_092121` and run `job.reproduce/reproduce_hep.ipynb`.


# Complete code

## Data curation: `job.curate_data`
Curate information for 74 diseases/traits: 
- Curate information for the 74 diseases: `job.curate_data/get_trait_list.ipynb`

Curate gene set (.gs) files:
- .gs file for 74 diseases: `job.curate_data/curate_gs_file.ipynb`
- .gs file for T cell signatures: `job.curate_data/curate_gs.tcell_signature.ipynb`
- .gs file for ploidy signatures: `job.curate_data/curate_ploidy_gs.ipynb`
- .gs file for zonation signatures: `job.curate_data/curate_zonation_gs.ipynb`
- .gs file for metabolic pathways: `job.curate_data/curate_gs.metabolic.ipynb`

Curate scRNA-seq data sets:
- TS FACS: 'job.curate_data/curate_ts_data.ipynb'
- Cano-Gamez & Soskic et al.: 'job.curate_data/curate_canogamez_tcell_data.ipynb'
- Nathan et al.: 'job.curate_data/curate_nathan_tcell_data.ipynb'
- Aizarani et al.: 'job.curate_data/curate_aizarani_liver_atlas_data.ipynb'
- Halpern & Shenhav et al.: 'job.curate_data/curate_halpern_mouse_liver_data.ipynb'
- Richter & Deligiannis et al.: 'job.curate_data/curate_richter_hepatocyte_data.ipynb'


## Compute scDRS scores: `job.compute_score`
- TMS FACS + 74 diseases: 'job.compute_score/compute_score.tms_facs_with_cov.magma_10kb_1000.sh'
- TMS FACS + T cell signatures: 'job.compute_score/compute_score.tms_facs_with_cov.tcell_sig.sh'
- TMS FACS + metabolic: 'job.compute_score/compute_score.tms_facs_with_cov.hep_metabolic.sh'
- TMS droplet + 74 diseases: 'job.compute_score/compute_score.tms_droplet_with_cov.magma_10kb_1000.sh'
- TS FACS + 74 diseases: 'job.compute_score/compute_score.ts_facs_with_cov.magma_10kb_1000.sh'


## Scehma (Fig. 1): `job.schema`
Make schematic figures.

## Simulation (Fig. 2): `job.simulation`
Data generation:
- Generate the TMS FACS 10K subsampled data and null gene sets: `job.simulation/generate_null_simulation_data.ipynb`
- Generate causal gene sets and perturbation configurations: `job.simulation/generate_causal_simulation_data.ipynb`

Compute results: 
- Compute scDRS scores for null simulations: `job.simulation/compute_simu_score.sh`
- Compute Seurat scores for null simulations: `job.simulation/compute_simu_score_scanpy.sh`
- Compute Vision scores for null simulations: `job.simulation/compute_simu_score_vision.sh`
- Compute VAM scores for null simulations: `job.simulation/call_R_vam.sh`
- Compute scores (scDRS/Seurat/Vision) for causal simulations (500 random causal cells): `job.simulation/compute_perturb_simu_score.sh`
- Compute scores (scDRS/Seurat/Vision) for causal simulations (B cells causal): `job.simulation/compute_perturb_simu_score_Bcell.sh`

Make figures:
- Make figures for null simulations: `job.simulation/make_figure.null_simulation.ipynb`
- Make figures for causal simulations (500 random causal cells): `job.simulation/make_figure.causal_simulation.ipynb`
- Make figures for causal simulations (B cells causal): `job.simulation/make_figure.causal_simulation_Bcell.ipynb`


## Cell type-level results (Fig. 3): `job.celltype_association`
XXX

## T cell example (Fig. 4): `job.case_tcell`
- Reprocess TMS T cells and assign effectorness gradients: `job.case_tcell/s1_reprocess_tms_tcell.ipynb`
- Main analysis: `job.case_tcell/s3_analysis_tcell.ipynb`
- Replication in Cano-Gamez & Soskic et al. and Nathan et al. data: `job.case_tcell/s4_analysis_tcell.replication.ipynb`
- Cluster-level LDSC-SEG analysis: `job.case_tcell/s5_compare_ldsc_cluster_4res.ipynb`
- Disease gene prioritization: `job.case_tcell/s6_gene_prioritization.ipynb`

## Neuron example (Fig. 5AB):  `job.ca1_pyramidal`
See details in `job.ca1_pyramidal/README.md`

## Hepatocyte example (Fig. 5CD): `job.case_hepatocyte`
- Reprocess TMS hepatocytes: `job.case_hepatocyte/s1_reprocess_tms_hep.ipynb`
- Main analysis: `job.case_hepatocyte/s3_analysis_hep.ipynb`
