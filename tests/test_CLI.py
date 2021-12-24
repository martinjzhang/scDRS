import scdrs
import os
import subprocess
import pandas as pd
import numpy as np
import tempfile
from .test_method_score_cell_main import compare_score_file


def test_score_cell_cli():
    """
    Test CLI `scdrs compute-score`
    """
    # Load toy data
    ROOT_DIR = scdrs.__path__[0]
    H5AD_FILE = os.path.join(ROOT_DIR, "data/toydata_mouse.h5ad")
    COV_FILE = os.path.join(ROOT_DIR, "data/toydata_mouse.cov")
    assert os.path.exists(H5AD_FILE), "built-in data toydata_mouse.h5ad missing"
    assert os.path.exists(COV_FILE), "built-in data toydata_mouse.cov missing"

    tmp_dir = tempfile.TemporaryDirectory()
    tmp_dir_path = tmp_dir.name
    dict_df_score = {}
    for gs_species in ["human", "mouse"]:
        gs_file = os.path.join(ROOT_DIR, f"data/toydata_{gs_species}.gs")
        # call compute_score.py
        cmds = [
            f"scdrs compute-score",
            f"--h5ad_file {H5AD_FILE}",
            "--h5ad_species mouse",
            f"--gs_file {gs_file}",
            f"--gs_species {gs_species}",
            f"--cov_file {COV_FILE}",
            "--ctrl_match_opt mean_var",
            "--n_ctrl 20",
            "--flag_filter_data False",
            "--weight_opt vs",
            "--flag_raw_count False",
            "--flag_return_ctrl_raw_score False",
            "--flag_return_ctrl_norm_score False",
            f"--out_folder {tmp_dir_path}",
        ]
        subprocess.check_call(" ".join(cmds), shell=True)
        dict_df_score[gs_species] = pd.read_csv(
            os.path.join(tmp_dir_path, f"toydata_gs_{gs_species}.score.gz"),
            sep="\t",
            index_col=0,
        )
    # consistency between human and mouse
    assert np.all(dict_df_score["mouse"].pval == dict_df_score["human"].pval)

    df_res = dict_df_score["mouse"]

    REF_COV_FILE = os.path.join(
        ROOT_DIR, "data/toydata_gs_mouse.ref_Ctrl20_CovConstCovariate.score.gz"
    )
    df_ref_res = pd.read_csv(REF_COV_FILE, sep="\t", index_col=0)
    compare_score_file(df_res, df_ref_res)
    tmp_dir.cleanup()
    return


def test_munge_gs_cli():
    """
    Test CLI `scdrs munge-gs`
    """

    tmp_dir = tempfile.TemporaryDirectory()
    tmp_dir_path = tmp_dir.name

    # pval_file and zscore_file
    temp_df = pd.DataFrame(
        data={
            "HEIGHT": [0.02, np.nan, 0.4],
            "BMI": [0.8, 0.02, np.nan],
        }
    )
    temp_df.index = ["OR4F5", "DAZ1", "BPY2B"]
    temp_df.to_csv(os.path.join(tmp_dir_path, "pval_file.tsv"), sep="\t", index=True)
    temp_df = pd.DataFrame(
        data={
            "GENE": ["OR4F5", "DAZ1", "BPY2B"],
            "HEIGHT": [2.0537, np.nan, 0.25335],
            "BMI": [-0.84162, 2.0537, np.nan],
        }
    )
    temp_df.to_csv(os.path.join(tmp_dir_path, "zscore_file.tsv"), sep="\t", index=False)

    dict_df_score = {}
    for input_file in ["pval_file", "zscore_file"]:
        for selection in [
            "--n-max 1",
            "--n-min 1 --n-max 3 --fdr 0.05",
            "--n-min 1 --n-max 3 --fwer 0.05",
        ]:
            # Call scdrs munge-gs
            input_file_path = os.path.join(tmp_dir_path, "%s.tsv" % input_file)
            output_file_path = os.path.join(tmp_dir_path, f"outfile.gs")
            cmds = [
                "scdrs munge-gs",
                f"--{input_file} {input_file_path}",
                f"--out-file {output_file_path}",
                "--weight zscore",
                selection,
            ]
            subprocess.check_call(" ".join(cmds), shell=True)
            temp_df = pd.read_csv(
                os.path.join(tmp_dir_path, f"outfile.gs"),
                sep="\t",
                index_col=0,
            )

            # Check results
            print('Generated .gs file:')
            print(temp_df)
            err_msg = "input_file=%s, %s" % (input_file, selection)
            assert list(temp_df.index) == ["BMI", "HEIGHT"], err_msg
            assert temp_df.loc["BMI", "GENESET"] == "DAZ1:2.0537", err_msg
            assert temp_df.loc["HEIGHT", "GENESET"] == "OR4F5:2.0537", err_msg

    tmp_dir.cleanup()

    return


def test_downstream_cli():
    """
    Test CLI `scdrs perform-downstream`

    1. --group-analysis cell_type
    2. --corr-analysis causal_variable,non_causal_variable,covariate
    3. --gene-analysis
    """

    # Load toy data
    ROOT_DIR = scdrs.__path__[0]
    H5AD_FILE = os.path.join(ROOT_DIR, "data/toydata_mouse.h5ad")
    SCORE_FILE = os.path.join(ROOT_DIR, "data/@.full_score.gz")
    REF_RES_DIR = os.path.join(ROOT_DIR, "data/")

    tmp_dir = tempfile.TemporaryDirectory()
    tmp_dir_path = tmp_dir.name
    for task in [
        "--group-analysis cell_type",
        "--corr-analysis causal_variable,non_causal_variable,covariate",
        "--gene-analysis",
    ]:
        # Call scdrs downstream
        cmds = [
            f"scdrs perform-downstream",
            f"--h5ad_file {H5AD_FILE}",
            f"--score-file {SCORE_FILE}",
            task,
            "--flag-filter-data False",
            "--flag-raw-count False",
            "--knn-n-neighbors 15",
            "--knn-n-pcs 20",
            f"--out-folder {tmp_dir_path}",
        ]
        subprocess.check_call(" ".join(cmds), shell=True)

    # Check consistency between computed results and reference results
    for prefix in ["toydata_gs_mouse.ref_Ctrl20_CovConstCovariate"]:
        for suffix in ["scdrs_group.cell_type", "scdrs_gene", "scdrs_cell_corr"]:
            res_path = os.path.join(tmp_dir_path, f"{prefix}.{suffix}")
            ref_res_path = os.path.join(REF_RES_DIR, f"{prefix}.{suffix}")
            df_res = pd.read_csv(res_path, sep="\t", index_col=0)
            df_ref_res = pd.read_csv(ref_res_path, sep="\t", index_col=0)
            print(df_res)
            assert np.allclose(
                df_res.values, df_ref_res.values
            ), '%s, %s'%(prefix, suffix)

    tmp_dir.cleanup()
    return