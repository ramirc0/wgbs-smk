"""Convert methylation bedMethyl files to strand-specific bigWig files."""

import argparse
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import polars as pl


def setup_logging():
    """Log to stdout.

    Returns
    -------
    logging.Logger
        Module logger.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(__name__)


def process_bed_to_bigwig(bed_gz_path, chromsizes_path, output_bw_path, strand, logger):
    """Write the methylation level of one strand of a bedMethyl file as a bigWig.

    A strand without sites gets an empty file instead of a bigWig.

    Parameters
    ----------
    bed_gz_path : str
        Gzipped bedMethyl file.
    chromsizes_path : str
        Chromosome sizes file.
    output_bw_path : pathlib.Path
        BigWig to write.
    strand : {"+", "-"}
        Strand to keep.
    logger : logging.Logger
        Logger for progress and errors.

    Returns
    -------
    bool
        True if the bigWig was written.
    """
    try:
        logger.info(f"Processing {bed_gz_path} for strand {strand} -> {output_bw_path}")

        bed_schema = {
            "chr": pl.Utf8,  # Contig or chromosome name
            "start": pl.UInt32,  # Start position (0 offset)
            "end": pl.UInt32,  # End position (1 offset)
            "name": pl.Utf8,  # Name of item
            "score": pl.UInt16,  # Score from 0-1000 (capped methylation informative coverage)
            "strand": pl.Utf8,  # Strand: +, - or ?
            "thick_start": pl.UInt32,  # Start of where display should be thick
            "thick_end": pl.UInt32,  # End of where display should be thick
            "color": pl.Utf8,  # Colour value (RGB)
            "coverage": pl.UInt32,  # Methylation informative coverage
            "methylation_level": pl.Float32,  # Percentage of reads showing methylation
            "ref_genotype": pl.Utf8,  # Reference genotype
            "sample_genotype": pl.Utf8,  # Sample genotype
            "quality_score": pl.Float32,  # Quality score for genotype call
        }

        df = pl.read_csv(
            bed_gz_path,
            separator="\t",
            has_header=False,
            skip_rows=1,
            schema=bed_schema,
        ).with_columns(
            pl.col("methylation_level") / 100,
        )

        logger.info(f"bedMethyl head:\n{df.head()}")

        filtered_df = df.filter(pl.col("strand") == strand).select(
            "chr", "start", "end", "methylation_level"
        )

        logger.info(f"stranded bedMethyl head:\n{filtered_df.head()}")

        if filtered_df.height == 0:
            logger.warning(f"No {strand} strand data found in {bed_gz_path}")
            with open(output_bw_path, "w") as _:
                pass
            return True

        with (
            tempfile.NamedTemporaryFile(
                mode="w+", suffix=".bedgraph"
            ) as unsorted_bedgraph,
            tempfile.NamedTemporaryFile(
                mode="w+", suffix=".bedgraph"
            ) as sorted_bedgraph,
        ):
            filtered_df.write_csv(
                unsorted_bedgraph.name, include_header=False, separator="\t"
            )
            unsorted_bedgraph.flush()

            cmd = (
                "sort -k1,1 -k2,2n {unsorted_bedgraph} > {sorted_bedgraph}"
                " && bedGraphToBigWig {sorted_bedgraph} {chromsizes_path} {output_bw_path}"
            )
            formatted_cmd = cmd.format(
                unsorted_bedgraph=unsorted_bedgraph.name,
                sorted_bedgraph=sorted_bedgraph.name,
                chromsizes_path=chromsizes_path,
                output_bw_path=output_bw_path,
            )

            logger.info(f"Running: {formatted_cmd}")

            try:
                result = subprocess.run(
                    formatted_cmd,
                    capture_output=True,
                    shell=True,
                    text=True,
                    check=True,
                )
            except subprocess.CalledProcessError as e:
                logger.error(f"bedGraphToBigWig failed for {output_bw_path}: {e}")
                logger.error(f"stdout: {e.stdout}")
                logger.error(f"stderr: {e.stderr}")
                return False

            if result.stderr:
                logger.warning(f"bedGraphToBigWig stderr: {result.stderr}")

            logger.info(f"Successfully created {output_bw_path}")
            return True

    except Exception as e:  # noqa: BLE001
        # Log and count the failure. main() exits 1 after trying every bigWig.
        logger.error(f"Error processing {bed_gz_path} for strand {strand}: {e}")
        return False


def build_parser():
    """Build the extract_signal.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the extract_signal.py flags.
    """
    parser = argparse.ArgumentParser(
        description="Convert methylation BED files to strand-specific BigWig files"
    )
    parser.add_argument(
        "--bed-gz",
        nargs=3,
        required=True,
        help="Input BED.gz files in order: CHG, CHH, CpG",
    )
    parser.add_argument("--chromsizes", required=True, help="Chromosome sizes file")
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for BigWig files"
    )
    parser.add_argument(
        "--barcode", required=True, help="Sample barcode for naming output files"
    )
    return parser


def main(argv=None):
    """Write one bigWig per context and strand, exiting 1 if any fails.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    logger = setup_logging()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bed_files = {"chg": args.bed_gz[0], "chh": args.bed_gz[1], "cpg": args.bed_gz[2]}

    output_bigwigs = {
        "cpg_neg": output_dir / f"{args.barcode}_cpg_neg.bw",
        "cpg_pos": output_dir / f"{args.barcode}_cpg_pos.bw",
        "chg_pos": output_dir / f"{args.barcode}_chg_pos.bw",
        "chg_neg": output_dir / f"{args.barcode}_chg_neg.bw",
        "chh_pos": output_dir / f"{args.barcode}_chh_pos.bw",
        "chh_neg": output_dir / f"{args.barcode}_chh_neg.bw",
    }

    tasks = [
        (bed_files["cpg"], args.chromsizes, output_bigwigs["cpg_neg"], "-"),
        (bed_files["cpg"], args.chromsizes, output_bigwigs["cpg_pos"], "+"),
        (bed_files["chg"], args.chromsizes, output_bigwigs["chg_pos"], "+"),
        (bed_files["chg"], args.chromsizes, output_bigwigs["chg_neg"], "-"),
        (bed_files["chh"], args.chromsizes, output_bigwigs["chh_pos"], "+"),
        (bed_files["chh"], args.chromsizes, output_bigwigs["chh_neg"], "-"),
    ]

    logger.info(f"Processing {len(tasks)} BigWig files")

    success_count = 0
    for task in tasks:
        logger.info(f"Starting: {task[2]}")
        if process_bed_to_bigwig(*task, logger):
            success_count += 1
            logger.info(f"Completed: {task[2]}")
        else:
            logger.error(f"Failed: {task[2]}")

    logger.info(
        f"Processing complete: {success_count}/{len(tasks)} files processed successfully"
    )

    if success_count != len(tasks):
        logger.error("Some files failed to process")
        sys.exit(1)
    else:
        logger.info("All files processed successfully")


if __name__ == "__main__":
    main()
