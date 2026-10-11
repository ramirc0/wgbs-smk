"""Unit tests for contig_pools.py, gemBS's contig pooling."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow/scripts"))

from contig_pools import contig_pools

LIMIT = 25_000_000
# GRCh38 chr1-22, X and Y lengths.
MAIN = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
    145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
    101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
    50818468, 156040895, 57227415,
]  # fmt: skip
NAMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def named(pools, sizes):
    """Replace contig indices by names.

    Parameters
    ----------
    pools : dict of str to list of int
        Output of `contig_pools`.
    sizes : list of tuple of (str, int)
        Its input.

    Returns
    -------
    dict of str to list of str
        Contig names per pool.
    """
    return {pool: [sizes[i][0] for i in members] for pool, members in pools.items()}


def test_grch38_gets_one_job_per_chromosome_and_one_pool():
    """Check the GRCh38 no-alt layout: chr1-22, X, Y alone, the rest in Pool@1."""
    small = [("chrM", 16569)] + [(f"chrUn_{i}", 68_000) for i in range(170)]
    sizes = list(zip(NAMES, MAIN, strict=True)) + small
    pools = named(contig_pools(sizes, LIMIT), sizes)
    assert list(pools) == [*NAMES, "Pool@1"]
    assert pools["Pool@1"] == [name for name, _ in small]


def test_small_contigs_fill_the_emptiest_pool_longest_first():
    """Check gemBS's greedy packing into ceil(total / limit) pools."""
    sizes = [("big", 30), ("a", 9), ("b", 8), ("c", 7), ("d", 3), ("e", 2)]
    pools = named(contig_pools(sizes, 10), sizes)
    assert pools == {
        "big": ["big"],
        "Pool@1": ["a"],
        "Pool@2": ["b", "e"],
        "Pool@3": ["c", "d"],
    }


def test_decoys_spill_into_a_second_pool():
    """Check that 40 Mb of small contigs make two balanced pools, each contig once."""
    small = [(f"decoy{i}", 40_000_000 // 2000) for i in range(2000)]
    sizes = list(zip(NAMES, MAIN, strict=True)) + small
    pools = contig_pools(sizes, LIMIT)
    members = pools["Pool@1"] + pools["Pool@2"]
    assert sorted(members) == list(range(24, len(sizes)))
    assert abs(len(pools["Pool@1"]) - len(pools["Pool@2"])) <= 1
    assert "Pool@3" not in pools


def test_no_small_contigs_means_no_pool():
    """Check that a reference of long contigs only gets no Pool@ entry."""
    assert list(contig_pools([("a", 10), ("b", 12)], 10)) == ["a", "b"]
