"""Merge the CG, CHG and CHH bedMethyl files into one sorted, indexed methylC file."""

import argparse
import os
import subprocess
from collections import defaultdict
from pathlib import Path

import polars as pl


def read_bedmethyl(bed_gz_path):
    """Read a gzipped bedMethyl file, skipping its track line.

    Parameters
    ----------
    bed_gz_path : str
        Gzipped bedMethyl file.

    Returns
    -------
    polars.DataFrame
        One row per site. `methylation_level` is a fraction, not a percentage.
    """
    bed_schema = {
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

    df = pl.read_csv(
        bed_gz_path, separator="\t", has_header=False, skip_rows=1, schema=bed_schema
    ).with_columns(
        pl.col("methylation_level") / 100,
    )

    return df


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
        "--threads", default=1, type=int, help="bgzip and tabix threads."
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    return parser


def main(argv=None):
    """Write `<outdir>/<barcode>_methylc.bed.gz` and its tabix index.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    methyl_type_to_methyl_file = {"CHH": args.chh, "CHG": args.chg, "CG": args.cpg}

    methyl_type_to_methyl_df = defaultdict()
    for methyl_type, methyl_file in methyl_type_to_methyl_file.items():
        df = read_bedmethyl(methyl_file)
        print(f"Loaded bedmethyl file for type {methyl_type}")
        df = df.with_columns(pl.lit(methyl_type).alias("type"))
        df = df.select(
            "chr", "start", "end", "type", "methylation_level", "strand", "coverage"
        )
        methyl_type_to_methyl_df[methyl_type] = df

    print("Loaded all bedmethyl files")

    methylc_df = pl.concat(list(methyl_type_to_methyl_df.values()))
    print("Concatenated bedmethyl files into methylc file")

    sorted_methylc_df = methylc_df.sort(["chr", "start"], descending=[False, False])
    print("Sorted methylc file by chr and start")

    methylc_file = args.outdir / f"{args.barcode}_methylc.bed"
    sorted_methylc_df.write_csv(methylc_file, separator="\t", include_header=False)
    print(f"Wrote methylc file: {methylc_file}")

    subprocess.check_call(f"bgzip -f -@ {args.threads} {methylc_file!s}", shell=True)
    subprocess.check_call(
        f"tabix -@ {args.threads} -p bed {methylc_file!s}.gz", shell=True
    )
    print("bgzipped and indexed methylc file")

    print("Finished")


if __name__ == "__main__":
    main()
