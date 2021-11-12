import tempfile
from contextlib import contextmanager
import os
import subprocess
import urllib.request
import pandas as pd
import numpy as np

# Step 1: download data
# Step 2: run test
# Step 3: check results


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def test_consistency():
    # download data
    tmp = tempfile.TemporaryDirectory()
    tmp_dir = tmp.name
    with cd(tmp_dir):
        urllib.request.urlretrieve(
            "https://figshare.com/ndownloader/files/30853717",
            "data.zip",
        )
        subprocess.check_call(
            "unzip data.zip && mkdir -p data/ && mv single_cell_data/zeisel_2015/* data/"
            " && rm data.zip && rm -r single_cell_data",
            shell=True,
        )
        # sample only the first trait
        subprocess.check_call(
            "cat data/geneset.gs | head -n 2 > data/geneset_subset.gs",
            shell=True,
        )
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        cmds = [
            f"python {cur_dir}/../compute_score.py",
            "--h5ad_file data/expr.h5ad",
            "--h5ad_species mouse",
            "--gs_file data/geneset_subset.gs",
            "--gs_species mouse",
            "--cov_file data/cov.tsv",
            "--flag_filter True",
            "--flag_raw_count True",
            "--flag_return_ctrl_raw_score False",
            "--flag_return_ctrl_norm_score True",
        ]
        subprocess.check_call("mkdir -p dense/", shell=True)
        subprocess.check_call(
            " ".join(cmds + ["--flag_sparse False --out_folder dense"]), shell=True
        )
        subprocess.check_call("mkdir -p sparse/", shell=True)
        subprocess.check_call(
            " ".join(cmds + ["--flag_sparse True --out_folder sparse"]), shell=True
        )
        a = (
            pd.read_csv("dense/PASS_MDD_Howard2019.score.gz", sep="\t")
            .iloc[:, 1:]
            .values
        )
        b = (
            pd.read_csv("sparse/PASS_MDD_Howard2019.score.gz", sep="\t")
            .iloc[:, 1:]
            .values
        )
        # check that the relative error between dense and sparse scores are within 0.1%
        assert np.allclose(a, b, rtol=1e-3)
    tmp.cleanup()
