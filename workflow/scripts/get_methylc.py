"""Merge the CG, CHG and CHH bedMethyl files into one sorted, indexed methylC file."""

import argparse
import subprocess
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


def contig_slices(df):
    """Split a frame sorted by contig into one slice per contig.

    Parameters
    ----------
    df : polars.DataFrame
        Frame with a `chr` column, each contig's rows together.

    Yields
    ------
    tuple of (str, polars.DataFrame)
        Contig name and its rows, a slice that shares the frame's memory.
    """
    offset = 0
    runs = df["chr"].rle().struct.unnest().select("value", "len")
    for contig, length in runs.iter_rows():
        yield contig, df.slice(offset, length)
        offset += length


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
    parser.add_argument(
        "--threads", default=1, type=int, help="Threads for bgzip and tabix."
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    return parser


def main(argv=None):
    """Write `<outdir>/<barcode>_methylc.bed.gz` and its tabix index.

    mextr sorts each file by contig name, then position, and no position is in
    two contexts. So each contig's rows of the three files, sorted by start,
    give the sorted methylC file.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    methylc = args.outdir / f"{args.barcode}_methylc.bed.gz"
    pieces = {}
    for context, bed_gz in (("CHH", args.chh), ("CHG", args.chg), ("CG", args.cpg)):
        for contig, piece in contig_slices(read_methylc(bed_gz, context)):
            pieces.setdefault(contig, []).append(piece)
    with open(methylc, "wb") as out:
        bgzip = subprocess.Popen(
            ["bgzip", "-@", str(args.threads)], stdin=subprocess.PIPE, stdout=out
        )
        # Code point order, as LC_ALL=C sort orders the UTF-8 names.
        for contig in sorted(pieces):
            pl.concat(pieces[contig]).sort("start").write_csv(
                bgzip.stdin, separator="\t", include_header=False
            )
        bgzip.stdin.close()
        if bgzip.wait() != 0:
            raise SystemExit("bgzip failed")
    subprocess.run(["tabix", "-@", str(args.threads), "-p", "bed", methylc], check=True)


if __name__ == "__main__":
    main()
