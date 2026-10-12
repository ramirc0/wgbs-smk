"""Merge the CG, CHG and CHH bedMethyl files into one sorted, indexed methylC file."""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import polars as pl

# mextr's bedMethyl columns, of which methylC keeps six.
BED_SCHEMA = {
    "chr": pl.Utf8,
    "start": pl.UInt32,
    "end": pl.UInt32,
    "name": pl.Utf8,
    "score": pl.UInt16,
    "strand": pl.Utf8,
    "thick_start": pl.UInt32,
    "thick_end": pl.UInt32,
    "color": pl.Utf8,
    "coverage": pl.UInt32,
    "methylation_level": pl.Float32,
    "ref_genotype": pl.Utf8,
    "sample_genotype": pl.Utf8,
    "quality_score": pl.Float32,
}


def read_methylc(bed_gz, context):
    """Read a bedMethyl file as methylC rows.

    Parameters
    ----------
    bed_gz : str
        Gzipped bedMethyl file with a track line.
    context : str
        Value of the `type` column, such as ``CG``.

    Returns
    -------
    polars.DataFrame
        Columns `chr`, `start`, `end`, `type`, `methylation_level` (a fraction),
        `strand` and `coverage`, in file order.
    """
    df = pl.read_csv(
        bed_gz,
        separator="\t",
        has_header=False,
        skip_rows=1,
        schema=BED_SCHEMA,
        columns=[0, 1, 2, 5, 9, 10],
    )
    return df.rechunk().select(
        "chr",
        "start",
        "end",
        pl.lit(context).alias("type"),
        pl.col("methylation_level") / 100,
        "strand",
        "coverage",
    )


def build_parser():
    """Build the get_methylc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the get_methylc.py flags.
    """
    parser = argparse.ArgumentParser(description="Build a methylC file per sample.")
    parser.add_argument("--barcode", required=True, help="Sample ID, the file prefix.")
    parser.add_argument("--chh", required=True, help="CHH bedMethyl file.")
    parser.add_argument("--chg", required=True, help="CHG bedMethyl file.")
    parser.add_argument("--cpg", required=True, help="CpG bedMethyl file.")
    parser.add_argument("--threads", default=1, type=int, help="Threads for bgzip.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--tmpdir", default=None, help="Directory for the parts.")
    return parser


def main(argv=None):
    """Write `<outdir>/<barcode>_methylc.bed.gz` and its tabix index.

    mextr sorts each file by contig name, then position. So a merge of the three
    files, with sort -m, gives the fully sorted methylC file.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    methylc = args.outdir / f"{args.barcode}_methylc.bed.gz"
    with tempfile.TemporaryDirectory(dir=args.tmpdir) as tmp:
        parts = []
        for context, bed_gz in (("CHH", args.chh), ("CHG", args.chg), ("CG", args.cpg)):
            part = Path(tmp, f"{context}.bed")
            read_methylc(bed_gz, context).write_csv(
                part, separator="\t", include_header=False
            )
            parts.append(part)
        merge = subprocess.Popen(
            ["sort", "-m", "-k1,1", "-k2,2n", *parts],
            stdout=subprocess.PIPE,
            env={**os.environ, "LC_ALL": "C"},
        )
        with open(methylc, "wb") as out:
            subprocess.run(
                ["bgzip", "-@", str(args.threads)],
                stdin=merge.stdout,
                stdout=out,
                check=True,
            )
        if merge.wait() != 0:
            raise SystemExit("sort -m failed")
    subprocess.run(["tabix", "-p", "bed", methylc], check=True)


if __name__ == "__main__":
    main()
