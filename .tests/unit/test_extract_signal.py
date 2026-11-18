"""Unit tests for the bigWig writes of extract_signal.py."""

import gzip
import sys
from pathlib import Path

import pyBigWig
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow/scripts"))

import extract_signal

# bigtools reads a bedGraph of 200 MB or more in parallel. A long contig name
# reaches that size in 1 M rows.
CONTIG = "c" * 200
LAST_CONTIG = CONTIG + "_random_unplaced_contig"
SITE = '{0}\t{1}\t{2}\t"S"\t5\t-\t{1}\t{2}\t0,255,0\t5\t20\tCG\tCG\t30\n'


def write_bed_gz(path, sites):
    """Write a gzipped bedMethyl file of minus-strand sites.

    Parameters
    ----------
    path : pathlib.Path
        The file to write.
    sites : iterable of tuple
        A contig and a start for each site.
    """
    with gzip.open(path, "wt", compresslevel=1) as f:
        f.write('track name="S"\n')
        f.writelines(SITE.format(contig, i, i + 1) for contig, i in sites)


def test_last_contig_with_one_site(tmp_path):
    """Check a large bedGraph whose last contig holds a single site.

    bigtools' parallel reader finds contig starts by bisecting the file. It
    misses a last contig of one row that is longer than the row before it, then
    stops with "File is not sorted".
    """
    chromsizes = tmp_path / "chromsizes.tsv"
    chromsizes.write_text(f"{CONTIG}\t2000000\n{LAST_CONTIG}\t1000\n")
    beds = [tmp_path / f"S_{context}.bed.gz" for context in ("chg", "chh", "cpg")]
    write_bed_gz(beds[0], [(CONTIG, 0)])
    write_bed_gz(beds[1], [(CONTIG, 0)])
    sites = [(CONTIG, i) for i in range(1_000_000)] + [(LAST_CONTIG, 5)]
    write_bed_gz(beds[2], sites)

    extract_signal.main(
        [
            "--bed-gz",
            *map(str, beds),
            "--chromsizes",
            str(chromsizes),
            "--output-dir",
            str(tmp_path),
            "--barcode",
            "S",
            "--tmpdir",
            str(tmp_path),
            "--threads",
            "8",
        ]
    )

    bw = pyBigWig.open(str(tmp_path / "S_cpg_neg.bw"))
    assert bw.intervals(LAST_CONTIG) == ((5, 6, pytest.approx(0.2)),)
