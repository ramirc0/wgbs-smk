"""Write the methylation level of mextr's bedMethyl files as one bigWig per context and strand."""

import argparse
import subprocess
import tempfile
from pathlib import Path

import polars as pl
import pyBigWig

# mextr's bedMethyl columns, of which the bigWigs use five.
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
STRANDS = {"pos": "+", "neg": "-"}


def read_levels(bed_gz):
    """Read the sites of a bedMethyl file with their methylation as a fraction.

    Parameters
    ----------
    bed_gz : str
        Gzipped bedMethyl file with a track line.

    Returns
    -------
    polars.DataFrame
        Columns `chr`, `start`, `end`, `strand` and `methylation_level`, in
        file order.
    """
    return (
        pl.read_csv(
            bed_gz,
            separator="\t",
            has_header=False,
            skip_rows=1,
            schema=BED_SCHEMA,
            columns=[0, 1, 2, 5, 10],
        )
        .rechunk()
        .with_columns(pl.col("methylation_level") / 100)
    )


def write_empty_bigwig(path, chromsizes):
    """Write a bigWig with the reference's contigs and no data.

    Parameters
    ----------
    path : pathlib.Path
        The bigWig to write.
    chromsizes : str
        Chromosome sizes file.
    """
    with open(chromsizes) as f:
        sizes = [(name, int(size)) for name, size in (line.split()[:2] for line in f)]
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(sizes)
    bw.close()


def build_parser():
    """Build the extract_signal.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the extract_signal.py flags.
    """
    parser = argparse.ArgumentParser(
        description="Write methylation bigWigs per context and strand."
    )
    parser.add_argument(
        "--bed-gz",
        nargs=3,
        required=True,
        help="mextr bedMethyl files in order: CHG, CHH, CpG.",
    )
    parser.add_argument("--chromsizes", required=True, help="Chromosome sizes file.")
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="Output directory."
    )
    parser.add_argument("--barcode", required=True, help="Sample ID, the file prefix.")
    parser.add_argument("--tmpdir", default=None, help="Directory for the bedGraphs.")
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for all bigWig writes; a quarter goes to each.",
    )
    return parser


def main(argv=None):
    """Write `<barcode>_<context>_<strand>.bw` for CpG, CHG and CHH on both strands.

    mextr writes its contigs in name order and each contig's sites in position
    order, which bigtools needs. A strand without sites gets a bigWig with no data.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.tmpdir) as tmp:
        jobs = []
        for context, bed_gz in zip(("chg", "chh", "cpg"), args.bed_gz, strict=True):
            levels = read_levels(bed_gz)
            for strand, sign in STRANDS.items():
                bigwig = args.output_dir / f"{args.barcode}_{context}_{strand}.bw"
                sites = levels.filter(pl.col("strand") == sign).drop("strand")
                if sites.is_empty():
                    write_empty_bigwig(bigwig, args.chromsizes)
                    continue
                bedgraph = Path(tmp, f"{context}_{strand}.bedgraph")
                sites.write_csv(bedgraph, include_header=False, separator="\t")
                threads = str(max(1, args.threads // 4))
                cmd = ["bigtools", "bedgraphtobigwig", "--inmemory", "-t", threads]
                jobs.append(subprocess.Popen([*cmd, bedgraph, args.chromsizes, bigwig]))
            del levels  # free it before reading the next context
        failed = [job.args for job in jobs if job.wait() != 0]
    if failed:
        raise SystemExit(f"bigtools failed: {failed}")


if __name__ == "__main__":
    main()
