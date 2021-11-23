import scdrs
import os
import subprocess
import pandas as pd
import numpy as np
from .test_score_cell import compare_score_file


def test_score_cell():

    # Load toy data
    ROOT_DIR = scdrs.__path__[0]
    H5AD_FILE = os.path.join(ROOT_DIR, "data/toydata_mouse.h5ad")
    COV_FILE = os.path.join(ROOT_DIR, "data/toydata_mouse.cov")
    assert os.path.exists(H5AD_FILE), "built-in data toydata_mouse.h5ad missing"
    assert os.path.exists(COV_FILE), "built-in data toydata_mouse.cov missing"

    import tempfile

    tmp_dir = tempfile.TemporaryDirectory()
    tmp_dir_path = tmp_dir.name
    dict_df_score = {}
    for gs_species in ["human", "mouse"]:
        gs_file = os.path.join(ROOT_DIR, f"data/toydata_{gs_species}.gs")
        # call compute_score.py
        cmds = [
            f"python {ROOT_DIR}/../compute_score.py",
            f"--h5ad_file {H5AD_FILE}",
            "--h5ad_species mouse",
            f"--gs_file {gs_file}",
            f"--gs_species {gs_species}",
            f"--cov_file {COV_FILE}",
            "--ctrl_match_opt mean_var",
            "--n_ctrl 20",
            "--flag_filter False",
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