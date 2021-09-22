# Experiments

## Data curation: `job.curate_data`
Curate disease information: 
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
See `./job.compute_score/readme.md` for details.


## Scehma (Fig. 1): `job.schema`

## Simulation (Fig. 2): `job.simulation`
Data generation:
- Generate the TMS FACS 10K subsampled data and null gene sets: `job.simulation/generate_null_simulation_data.ipynb`
- Generate causal gene sets and perturbation configurations: `job.simulation/generate_causal_simulation_data.ipynb`

Compute results: 
- Compute scDRS scores for null simulation: `job.simulation/compute_simu_score.sh`
- Compute Seurat scores for null simulation: `job.simulation/compute_simu_score_scanpy.sh`
- Compute Vision scores for null simulation: `job.simulation/compute_simu_score_vision.sh`
- Compute VAM scores for null simulation: `job.simulation/call_R_vam.sh`
- Compute scores (scDRS/Seurat/Vision) for causal simulation (500 random causal cells): `job.simulation/compute_perturb_simu_score.sh`
- Compute scores (scDRS/Seurat/Vision) for causal simulation (B cell causal): `job.simulation/compute_perturb_simu_score_Bcell.sh`

Make figures:
- Make figures for null simulations: `job.simulation/make_figure.null_simulation.ipynb`
- Make figures for causal simulations with 500 random causal cells: `job.simulation/make_figure.causal_simulation.ipynb`
- Make figures for causal simulations with 528 B cells as causal cells: `job.simulation/make_figure.causal_simulation_Bcell.ipynb`


## Cell type-level results (Fig. 3): `job.celltype_association`

## T cell example (Fig. 4): `job.case_tcell`
- 

## Neuron example (Fig. 5AB):  `job.ca1_pyramidal`
## Hepatocyte example (Fig. 5CD): `job.case_hepatocyte`
