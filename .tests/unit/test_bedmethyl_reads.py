"""Unit tests for the bedMethyl readers of get_methylc.py and extract_signal.py."""

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow/scripts"))

from extract_signal import read_levels
from get_methylc import read_methylc

SITE = 'chr1\t{0}\t{1}\t"S"\t5\t+\t{0}\t{1}\t0,255,0\t5\t20\tCTT\tCTT\t30\n'


@pytest.fixture
def bed_gz(tmp_path):
    """Write a bedMethyl file big enough for polars to read it in several chunks.

    Returns
    -------
    pathlib.Path
        Gzipped bedMethyl file with a track line.
    """
    path = tmp_path / "S_chh.bed.gz"
    with gzip.open(path, "wt") as f:
        f.write('track name="S"\n')
        f.writelines(SITE.format(i, i + 1) for i in range(100_000))
    return path


@pytest.mark.parametrize(
    "read", [lambda p: read_methylc(p, "CHH"), read_levels], ids=["methylc", "levels"]
)
def test_reads_are_one_chunk(bed_gz, read):
    """Check that a read gives one chunk per column.

    polars writes a CSV of a frame in 140,000 chunks 25 times slower than one of
    the same frame in one chunk.
    """
    assert read(bed_gz).n_chunks() == 1
