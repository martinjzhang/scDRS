import pytest
import numpy as np
import scipy as sp

import scdrs


def test_select_ctrl_geneset():
    """
    Test scdrs.method._select_ctrl_geneset
    """
    return


def test_compute_raw_score():
    """
    Test scdrs.method._compute_raw_score
    """
    return


def test_correct_background():
    """
    Test scdrs.method._correct_background
    """
    return


def test_get_p_from_empi_null():
    """
    Test scdrs.method._get_p_from_empi_null
    """
    v_t = [0, 1]
    v_t_null = [0.5, 0.6, 0.7]
    v_p = scdrs.method._get_p_from_empi_null(v_t, v_t_null)
    assert np.absolute(v_p[0] - 1) < 0.001, "First MC-p should be 1"
    assert np.absolute(v_p[1] - 0.25) < 0.001, "Second MC-p should be 1/(3+1)=0.25"
    return


def test_get_rank():
    """
    Test scdrs.method._get_rank
    """
    return
