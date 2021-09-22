# Experiments

## Data curation: `job.curate_data`
Curate scRNA-seq data sets:

Curate other files: 
- Curate information for the 74 diseases: `job.curate_data/get_trait_list.ipynb`

## Compute scDRS scores: `job.compute_score`
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
