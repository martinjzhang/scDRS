import pytest
import numpy as np
import scipy as sp
import pandas as pd
import scdrs
import tempfile


def test_gearys_c():
    """
    gearys_c
    """

    return


def test_compute_gearysc_significance():
    """
    compute_gearysc_significance
    """
    return


def test_load_homolog_mapping():
    """
    Test scdrs.util.load_homolog_mapping
    """
    dict_mouse_to_human = scdrs.util.load_homolog_mapping(
        src_species="mmusculus", dst_species="hsapiens"
    )
    # check only take the first gene
    g_mouse, g_human = list(dict_mouse_to_human.items())[0]
    assert g_human.isupper(), "Human gene should be upper case '%s'" % g_huamn
    assert g_mouse[0].isupper() & g_mouse[1:].islower(), (
        "Mouse gene should be lower case '%s'" % g_mouse
    )

    dict_human_to_mouse = scdrs.util.load_homolog_mapping(
        src_species="hsapiens", dst_species="mmusculus"
    )
    # check only take the first gene
    g_human, g_mouse = list(dict_human_to_mouse.items())[0]
    assert g_human.isupper(), "Human gene should be upper case '%s'" % g_huamn
    assert g_mouse[0].isupper() & g_mouse[1:].islower(), (
        "Mouse gene should be lower case '%s'" % g_mouse
    )
    return


def test_load_gs():
    """
    Test scdrs.util.load_gs
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Uniform weights + same species
        df_gs1 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["ACADM,AGA"]})
        df_gs1.to_csv(tmp_dir + "/test.gs", index=False, sep="\t")
        dict_gs1 = scdrs.util.load_gs(tmp_dir + "/test.gs")
        assert dict_gs1["TRAIT1"] == (
            ["ACADM", "AGA"],
            [1.0, 1.0],
        ), "Uniform weights + same species"
        
         # Uniform weights + same species + duplicated values
        df_gs1 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["ACADM,ACADM,AGA"]})
        df_gs1.to_csv(tmp_dir + "/test.gs", index=False, sep="\t")
        dict_gs1 = scdrs.util.load_gs(tmp_dir + "/test.gs")
        assert dict_gs1["TRAIT1"] == (
            ["ACADM", "AGA"],
            [1.0, 1.0],
        ), "Uniform weights + same species + duplicated values"

        # Continuous weights + same species
        df_gs2 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["Acadm:-1.4,Aga:4.5"]})
        df_gs2.to_csv(tmp_dir + "/test.gs", index=False, sep="\t")
        dict_gs1 = scdrs.util.load_gs(tmp_dir + "/test.gs")
        assert dict_gs1["TRAIT1"] == (
            ["Acadm", "Aga"],
            [-1.4, 4.5],
        ), "Continuous weights + same species"
        
        # Continuous weights + same species + duplicated values (the later one was read)
        df_gs2 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["Acadm:-1.4,Acadm:2,Aga:4.5"]})
        df_gs2.to_csv(tmp_dir + "/test.gs", index=False, sep="\t")
        dict_gs1 = scdrs.util.load_gs(tmp_dir + "/test.gs")
        assert dict_gs1["TRAIT1"] == (
            ["Acadm", "Aga"],
            [2.0, 4.5],
        ), "Continuous weights + same species + duplicated values"

        # Continuous weights + human to mouse
        df_gs2 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["ACADM:2.30,AGA:4.5"]})
        df_gs2.to_csv(tmp_dir + "/test.gs", index=False, sep="\t")
        dict_gs1 = scdrs.util.load_gs(
            tmp_dir + "/test.gs", src_species="human", dst_species="mouse"
        )
        assert dict_gs1["TRAIT1"] == (
            ["Acadm", "Aga"],
            [2.3, 4.5],
        ), "Continuous weights + human to mouse"

        # Continuous weights + mouse to human
        df_gs2 = pd.DataFrame({"TRAIT": ["TRAIT1"], "GS": ["Acadm:2.3,Aga:4.5"]})
        df_gs2.to_csv(tmp_dir + "/test.gs", index=False, sep="\t")
        dict_gs1 = scdrs.util.load_gs(
            tmp_dir + "/test.gs", src_species="mouse", dst_species="human"
        )
        assert dict_gs1["TRAIT1"] == (
            ["ACADM", "AGA"],
            [2.3, 4.5],
        ), "Continuous weights + mouse to human"
    return
