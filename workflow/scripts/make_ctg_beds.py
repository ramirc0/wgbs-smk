"""Write one contig BED per main chromosome and one for all other contigs."""

import argparse
import os

import polars as pl


def parse_contig_sizes(contig_file: str):
    """Read a contig sizes file as whole-contig intervals.

    Parameters
    ----------
    contig_file : str
        Two-column contig sizes file.

    Returns
    -------
    polars.DataFrame
        Columns `chrom`, `start` (always 0) and `end`.
    """
    df = pl.read_csv(
        contig_file, separator="\t", has_header=False, new_columns=["chrom", "end"]
    )
    df = df.with_columns(pl.lit(0).alias("start")).select("chrom", "start", "end")
    return df


def create_bed_files(df: pl.DataFrame, output_dir: str, barcode: str):
    """Write `<barcode>_<chrom>_ctgs.bed` per main chromosome and `_Pool_ctgs.bed`.

    Main chromosomes are chr1-chr22, chrX and chrY. Absent ones get no file.
    The Pool file holds every other contig.

    Parameters
    ----------
    df : polars.DataFrame
        Contig intervals from `parse_contig_sizes`.
    output_dir : str
        Output directory.
    barcode : str
        Sample ID, the file prefix.
    """
    os.makedirs(output_dir, exist_ok=True)
    main_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    for chrom in main_chroms:
        chrom_data = df.filter(pl.col("chrom") == chrom)
        if not chrom_data.is_empty():
            bed_file = os.path.join(output_dir, f"{barcode}_{chrom}_ctgs.bed")
            chrom_data.write_csv(bed_file, include_header=False, separator="\t")
            print(f"Created: {bed_file}")

    alt_contigs = df.filter(~(pl.col("chrom").is_in(main_chroms)))
    if not alt_contigs.is_empty():
        pool_file = os.path.join(output_dir, f"{barcode}_Pool_ctgs.bed")
        alt_contigs.write_csv(pool_file, include_header=False, separator="\t")
        print(f"Created: {pool_file}")
        print(f"Pool contains {len(alt_contigs)} alternative contigs")


def build_parser():
    """Build the make_ctg_beds.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the make_ctg_beds.py flags.
    """
    parser = argparse.ArgumentParser(
        description="Generate contig BED files from contig sizes"
    )
    parser.add_argument("--chromsizes", help="Path to contig sizes file")
    parser.add_argument("--output_dir", help="Output directory for BED files")
    parser.add_argument("--barcode", help="Sample barcode/name")
    return parser


def main(argv=None):
    """Read the contig sizes and write the contig BEDs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    print(f"Reading chromsizes from: {args.chromsizes}")
    df = parse_contig_sizes(args.chromsizes)
    print(f"Found {df.shape[0]} contigs")

    create_bed_files(df, args.output_dir, args.barcode)
    print("Done!")


if __name__ == "__main__":
    main()
