import pytest
import numpy as np
import scipy as sp
import pandas as pd
import scdrs.util as util


"""
    Example
"""


def test_example():
    x = 1
    assert x == 1
    return


"""
    gearys_c
"""


def test_gearys_c():
    return


"""
    compute_gearysc_significance
"""


def test_compute_gearysc_significance():
    return


def test_gs_util():

    # Test load_mouse_human_homologs
    dict_mouse_to_human = util.load_mouse_human_homologs(
        src_species="mmusculus", dst_species="hsapiens"
    )
    # check only take the first gene
    g_mouse, g_human = list(dict_mouse_to_human.items())[0]
    assert g_human.isupper()
    assert g_mouse[0].isupper() & g_mouse[1:].islower()

    dict_human_to_mouse = util.load_mouse_human_homologs(
        src_species="hsapiens", dst_species="mmusculus"
    )
    # check only take the first gene
    g_human, g_mouse = list(dict_human_to_mouse.items())[0]
    assert g_human.isupper()
    assert g_mouse[0].isupper() & g_mouse[1:].islower()

    # Test load_gs
    import tempfile

    with tempfile.TemporaryDirectory() as tmp_dir:
        # check uniform weights
        df_gs1 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["GENE1,GENE2"]})
        df_gs1.to_csv(tmp_dir + "/gs1.csv", index=False, sep="\t")
        dict_gs1 = util.load_gs(tmp_dir + "/gs1.csv")
        assert dict_gs1["TRAIT1"][0] == ["GENE1", "GENE2"]
        assert dict_gs1["TRAIT1"][1] == [1.0, 1.0]

        # check quantitative weights
        df_gs2 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["GENE1:2.3,GENE2:4.5"]})
        df_gs2.to_csv(tmp_dir + "/gs2.csv", index=False, sep="\t")
        dict_gs1 = util.load_gs(tmp_dir + "/gs2.csv")
        assert dict_gs1["TRAIT1"][0] == ["GENE1", "GENE2"]
        assert dict_gs1["TRAIT1"][1] == [2.3, 4.5]