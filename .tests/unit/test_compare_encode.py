"""Unit tests for the validation's compare_encode.py."""

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))

from compare_encode import read_cpg


def test_read_cpg_keeps_five_columns(tmp_path):
    """Check that a mextr bedMethyl with a track line reads as five columns."""
    bed = tmp_path / "S_cpg.bed.gz"
    with gzip.open(bed, "wt") as f:
        f.write('track name = "S" description="S" visibility=2 itemRgb="On"\n')
        f.write('chr1\t32\t33\t"S"\t1\t-\t32\t33\t255,0,0\t7\t42.86\tCG\tCG\t2\n')
    df = read_cpg(bed)
    assert df.columns == ["chrom", "start", "strand", "coverage", "percent"]
    assert df.row(0) == ("chr1", 32, "-", 7, df["percent"][0])
    assert abs(df["percent"][0] - 42.86) < 1e-4
