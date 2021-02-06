"""Test covid-spike-classification core functions."""

import pytest

from covid_spike_classification import core


def test_parse_pileup_match():
    pileup = """NC_045512	22919	T	1	.	G
NC_045512	22920	A	1	.	E
NC_045512	22921	T	1	.	O
"""
    expected = ("TAT", "TAT", [38, 36, 46])
    assert expected == core.parse_pileup(pileup)


def test_parse_pileup_mismatch():
    pileup = """NC_045512	22919	T	1	.	G
NC_045512	22920	A	1	G	E
NC_045512	22921	T	1	.	O
"""
    expected = ("TAT", "TGT", [38, 36, 46])
    assert expected == core.parse_pileup(pileup)


def test_parse_pileup_deletion_end():
    pileup = """NC_045512	22919	T	1	.	G
NC_045512	22920	A	1	.	E
NC_045512	22921	T	1	.*G	O
"""
    expected = ("TAT", "TAT", [38, 36, 46])
    assert expected == core.parse_pileup(pileup)


def test_parse_pileup_read_start():
    pileup = """NC_045512	22919	T	1	^M.	G
NC_045512	22920	A	1	.	E
NC_045512	22921	T	1	.	O
"""
    expected = ("TAT", "TAT", [38, 36, 46])
    assert expected == core.parse_pileup(pileup)


def test_parse_pileup_too_short():
    pileup = """NC_045512	22919	T	1	.	G
NC_045512	22920	A	1	.	E
"""
    with pytest.raises(core.PileupFailedError):
        _ = core.parse_pileup(pileup)
