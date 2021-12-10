import scanpy as sc
from anndata import read_h5ad
import pandas as pd
import numpy as np
import scipy as sp
import os
import fnmatch
import time
import argparse
from statsmodels.stats.multitest import multipletests

# Inhouse tools
import scdrs.util as util
import scdrs.data_loader as dl
import scdrs.method as md


"""
# Fixit


# Todo
- Implement a memory efficient version
- "gene_weight" argument needs to be tested 

# Finished
- Add --n_ctrl (default value 500) 
- Add --cov_file option to regress out covariates stored in COV_FILE before feeding into the score function 
- Add --ctrl_match_opt='mean_var': use mean- and var- matched control genes 
- Change name from scTRS to scdrs (072721)
- Fixed: Warning for compute_score: Trying to set attribute `.X` of view, copying. (did: v_norm_score = v_raw_score.copy())

"""

VERSION = "0.0.1"
VERSION = "beta"


def main(args):
    sys_start_time = time.time()

    MASTHEAD = "******************************************************************************\n"
    MASTHEAD += "* scDRS downsteam analyses \n"
    MASTHEAD += "* Version %s\n" % VERSION
    MASTHEAD += "* Martin Jinye Zhang and Kangcheng Hou\n"
    MASTHEAD += "* HSPH / Broad Institute / UCLA\n"
    MASTHEAD += "* MIT License\n"
    MASTHEAD += "******************************************************************************\n"

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################
    H5AD_FILE = args.h5ad_file
    SCORE_FILE = args.score_file
    CELLTYPE_LIST = [] if args.cell_type is None else args.cell_type.split(",")
    VARIABLE_LIST = [] if args.cell_variable is None else args.cell_variable.split(",")
    FLAG_GENE = args.flag_gene == "True"
    FLAG_FILTER = args.flag_filter == "True"
    FLAG_RAW_COUNT = args.flag_raw_count == "True"
    OUT_FOLDER = args.out_folder

    header = MASTHEAD
    header += "Call: ./compute_downstream.py \\\n"
    header += "--h5ad_file %s\\\n" % H5AD_FILE
    header += "--score_file %s\\\n" % SCORE_FILE
    header += "--cell_type %s\\\n" % args.cell_type
    header += "--cell_variable %s\\\n" % args.cell_variable
    header += "--flag_gene %s\\\n" % FLAG_GENE
    header += "--flag_filter %s\\\n" % FLAG_FILTER
    header += "--flag_raw_count %s\\\n" % FLAG_RAW_COUNT
    header += "--out_folder %s\n" % OUT_FOLDER
    print(header)

    ###########################################################################################
    ######                                     Load data                                 ######
    ###########################################################################################
    print("Load data:")

    # Load .h5ad file
    adata = read_h5ad(H5AD_FILE)
    if FLAG_FILTER:
        sc.pp.filter_cells(adata, min_genes=250)
        sc.pp.filter_genes(adata, min_cells=50)
    if FLAG_RAW_COUNT:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
    print(
        "--h5ad_file loaded: n_cell=%d, n_gene=%d (sys_time=%0.1fs)"
        % (adata.shape[0], adata.shape[1], time.time() - sys_start_time)
    )

    # Check CELLTYPE_LIST and VARIABLE_LIST
    temp_list = [x for x in CELLTYPE_LIST + VARIABLE_LIST if x not in adata.obs.columns]
    if len(temp_list) > 0:
        raise ValueError(
            "Following columns not in adata.obs.columns: %s" % ",".join(temp_list)
        )
    else:
        print("cell_type and cell_variable are in adata.obs.columns")

    # Load score file
    score_file_pattern = SCORE_FILE.split(os.path.sep)[-1]
    score_dir = SCORE_FILE.replace(os.path.sep + score_file_pattern, "")
    score_file_list = [
        x
        for x in os.listdir(score_dir)
        if fnmatch.fnmatch(x, score_file_pattern.replace("@", "*"))
    ]
    print("Infer score_dir=%s" % score_dir)
    print("Find %s score_files: %s" % (len(score_file_list), ",".join(score_file_list)))
    dic_score = {}
    for score_file in score_file_list:
        temp_df = pd.read_csv(
            score_dir + os.path.sep + score_file, sep="\t", index_col=0
        )
        n_cell_overlap = len(set(adata.obs_names) & set(temp_df.index))
        if n_cell_overlap < 0.1 * adata.shape[0]:
            print(
                "WARNING: %s skipped, %d/%d cells in adata"
                % (score_file, n_cell_overlap, adata.shape[0])
            )
        else:
            dic_score[score_file.replace(".full_score.gz", "")] = temp_df.copy()

    print(
        "--score_file loaded: n_trait=%d, (sys_time=%0.1fs)"
        % (len(dic_score), time.time() - sys_start_time)
    )
    print("")

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################
    STR_ANALYSIS = "Perform downstream analyses:"
    i = 1
    for ct in CELLTYPE_LIST:
        STR_ANALYSIS += "\n%d. Cell type-level analysis using %s" % (i, ct)
        STR_ANALYSIS += ": results in @.scdrs_ct.%s" % (ct)
        i += 1
    if len(VARIABLE_LIST) > 0:
        STR_ANALYSIS += "\n%d. Variable-disease correlation analysis for (%s)" % (
            i,
            ",".join(VARIABLE_LIST),
        )
        STR_ANALYSIS += ": results in @.scdrs_var"
        i += 1
    if FLAG_GENE is True:
        STR_ANALYSIS += "\n%d. Disease gene prioritization" % i
        STR_ANALYSIS += ": results in @.scdrs_gene"
    print(STR_ANALYSIS)

    # Compute connectivities if need to do cell type-level analysis
    if (len(CELLTYPE_LIST) > 0) & ("connectivities" not in adata.obsp):
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        print(
            "Compute connectivities with `sc.pp.neighbors` because `connectivities` is not found in adata.obsp"
        )

    # A separate file for each trait
    for trait in dic_score.keys():
        cell_list = sorted(set(adata.obs_names) & set(dic_score[trait].index))
        control_list = [
            x for x in dic_score[trait].columns if x.startswith("ctrl_norm_score")
        ]
        n_ctrl = len(control_list)
        df_reg = adata.obs.loc[cell_list, CELLTYPE_LIST + VARIABLE_LIST].copy()
        df_reg = df_reg.join(
            dic_score[trait].loc[cell_list, ["norm_score"] + control_list]
        )

        # Cell type-disease analysis: association+heterogeneity
        for ct_col in CELLTYPE_LIST:
            ct_list = sorted(set(adata.obs[ct_col]))
            col_list = [
                "n_cell",
                "n_ctrl",
                "assoc_mcp",
                "assoc_mcz",
                "hetero_mcp",
                "hetero_mcz",
            ]
            df_res = pd.DataFrame(index=ct_list, columns=col_list, dtype=np.float32)
            # Basic info
            for ct in ct_list:
                ct_cell_list = list(df_reg.index[df_reg[ct_col] == ct])
                df_res.loc[ct, ["n_cell", "n_ctrl"]] = [len(ct_cell_list), n_ctrl]
            # Association
            for ct in ct_list:
                ct_cell_list = list(df_reg.index[df_reg[ct_col] == ct])
                score_q95 = np.quantile(df_reg.loc[ct_cell_list, "norm_score"], 0.95)
                v_ctrl_score_q95 = np.quantile(
                    df_reg.loc[ct_cell_list, control_list], 0.95, axis=0
                )
                mc_p = ((v_ctrl_score_q95 >= score_q95).sum() + 1) / (
                    v_ctrl_score_q95.shape[0] + 1
                )
                mc_z = (score_q95 - v_ctrl_score_q95.mean()) / v_ctrl_score_q95.std()
                df_res.loc[ct, ["assoc_mcp", "assoc_mcz"]] = [mc_p, mc_z]
            # Heterogeneity
            # subset to common set of cells
            df_rls = md.test_gearysc(
                adata[cell_list], df_reg.loc[cell_list, :], groupby=ct_col
            )
            for ct in ct_list:
                mc_p, mc_z = df_rls.loc[ct, ["pval", "zsc"]]
                df_res.loc[ct, ["hetero_mcp", "hetero_mcz"]] = [mc_p, mc_z]

            df_res.to_csv(
                os.path.join(
                    OUT_FOLDER, "%s.scdrs_ct.%s" % (trait, ct_col.replace(" ", "_"))
                ),
                sep="\t",
                index=True,
            )
            print(
                "%s: cell type-level analysis with label=%s (sys_time=%0.1fs)"
                % (trait, ct_col, time.time() - sys_start_time)
            )

        # Variable-disease correlation
        if len(VARIABLE_LIST) > 0:
            col_list = ["n_ctrl", "corr_mcp", "corr_mcz"]
            df_res = pd.DataFrame(
                index=VARIABLE_LIST, columns=col_list, dtype=np.float32
            )
            for var_col in VARIABLE_LIST:
                corr_ = np.corrcoef(df_reg[var_col], df_reg["norm_score"])[0, 1]
                v_corr_ = [
                    np.corrcoef(df_reg[var_col], df_reg["ctrl_norm_score_%d" % x])[0, 1]
                    for x in np.arange(n_ctrl)
                ]
                v_corr_ = np.array(v_corr_)
                mc_p = ((v_corr_ >= corr_).sum() + 1) / (v_corr_.shape[0] + 1)
                mc_z = (corr_ - v_corr_.mean()) / v_corr_.std()
                df_res.loc[var_col] = [n_ctrl, mc_p, mc_z]
            df_res.to_csv(
                os.path.join(OUT_FOLDER, "%s.scdrs_var" % trait), sep="\t", index=True
            )
            print(
                "%s: cell-level variable-disease correlation analysis (sys_time=%0.1fs)"
                % (trait, time.time() - sys_start_time)
            )

        # Gene prioritization
        if FLAG_GENE is True:
            mat_expr = adata[df_reg.index].X.copy()
            v_corr = md._pearson_corr(mat_expr, df_reg["norm_score"].values)
            df_res = pd.DataFrame(
                index=adata.var_names, columns=["CORR", "RANK"], dtype=np.float32
            )
            df_res["CORR"] = v_corr
            df_res.sort_values("CORR", ascending=False, inplace=True)
            df_res["RANK"] = np.arange(df_res.shape[0])
            df_res.to_csv(
                os.path.join(OUT_FOLDER, "%s.scdrs_gene" % trait), sep="\t", index=True
            )
            print(
                "%s: disease gene prioritization (sys_time=%0.1fs)"
                % (trait, time.time() - sys_start_time)
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute score")

    parser.add_argument("--h5ad_file", type=str, required=True)
    parser.add_argument(
        "--score_file",
        type=str,
        required=True,
        help="@.full_score.gz where @ denotes trait names",
    )
    parser.add_argument(
        "--cell_type",
        type=str,
        required=False,
        default=None,
        help="Comma-seprated coloumn names for cell types/tissues, "
        "used for assessing cell type-disease association and "
        "within-cell type disease association heterogeneity",
    )
    parser.add_argument(
        "--cell_variable",
        type=str,
        required=False,
        default=None,
        help="Comma-seprated coloumn names for cell-level variables, "
        "used for associating cell-level variables to disease scores",
    )
    parser.add_argument(
        "--flag_gene",
        type=str,
        required=False,
        default=False,
        help="If True, perform gene prioritization",
    )
    parser.add_argument(
        "--flag_filter",
        type=str,
        required=False,
        default="True",
        help="If to apply cell and gene filters to the h5ad_file data",
    )
    parser.add_argument(
        "--flag_raw_count",
        type=str,
        required=False,
        default="True",
        help="If True, apply size factor normalization and log1p transformation",
    )
    parser.add_argument(
        "--out_folder",
        type=str,
        required=True,
        help="Save file at out_folder/trait.scdrs_res",
    )

    args = parser.parse_args()

    main(args)