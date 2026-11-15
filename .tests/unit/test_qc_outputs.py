"""Unit tests for qc_outputs.py, the readers behind the QC report."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow/scripts"))

from qc_outputs import (
    conversion_pct,
    depth_ratio,
    mean_methylation_pct,
    read_call_reports,
    sites_at_depth,
)


def call_report(path, profile, methylation, depths, passed):
    """Write a bs_call report with the fields the readers use.

    Parameters
    ----------
    path : pathlib.Path
        JSON file to write.
    profile : list of list of int
        ``NonCpGreadProfile`` rows.
    methylation : list of float
        ``PassedRefCpg`` histogram.
    depths : dict of str to int
        ``RefCpG`` coverage histogram.
    passed : int
        Passed reference CpGs.
    """
    stats = {
        "methylation": {"NonCpGreadProfile": profile, "PassedRefCpg": methylation},
        "coverage": {"RefCpG": depths},
        "RefCpG": {"All": passed + 1, "Passed": passed},
    }
    path.write_text(json.dumps({"totalStats": stats}))


def test_call_reports_add_up_across_pools(tmp_path):
    """Check that the pools' profiles, histograms and counts are summed."""
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    call_report(a, [[1, 2, 3, 4]], [1.0, 0.0], {"3": 5}, 7)
    call_report(b, [[1, 1, 1, 1], [9, 9, 9, 9]], [0.5, 2.0], {"3": 1, "12": 2}, 3)
    calls = read_call_reports([a, b])
    assert calls["profile"] == [[2, 3, 4, 5], [9, 9, 9, 9]]
    assert calls["methylation"] == [1.5, 2.0]
    assert calls["depths"] == {3: 6, 12: 2}
    assert calls["passed_cpgs"] == 10


def test_conversion_counts_converted_cs_at_non_cpg_sites():
    """Check that conversion is 1 minus gemBS's plotted non-conversion, a / (a + b)."""
    assert conversion_pct([[1, 99, 50, 0], [0, 100, 50, 1]]) == pytest.approx(99.5)
    assert conversion_pct([[0, 0, 0, 0]]) is None


def test_mean_methylation_reads_the_percent_bins():
    """Check that bin i of the 101 bins stands for i percent methylation."""
    histogram = [0.0] * 101
    histogram[0] = 1.0
    histogram[80] = 3.0
    assert mean_methylation_pct(histogram) == pytest.approx(60.0)


def test_sites_at_depth_counts_from_the_threshold_up():
    """Check that sites at exactly the threshold count."""
    assert sites_at_depth({3: 6, 10: 2, 40: 1}, 10) == 3


def test_depth_ratio_scales_by_contig_length():
    """Check a contig's reads per base against the autosomes'."""
    reads = {"chr1": 100, "chr2": 300, "chrX": 50, "chrM": 999}
    sizes = {"chr1": 10, "chr2": 30, "chrX": 10, "chrM": 1}
    assert depth_ratio(reads, sizes, "X") == pytest.approx(0.5)
    assert depth_ratio(reads, sizes, "Y") is None
    assert depth_ratio({"ctg1": 1, "ctgX": 1}, {"ctg1": 1, "ctgX": 1}, "X") is None
