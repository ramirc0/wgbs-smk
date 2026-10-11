"""Write the mapping QC JSON: the gem-mapper reports, samtools stats and average coverage.

Coverage = (mapped reads * average read length) / genome size. The genome size is
the sum of the chromosome sizes. Adapted from
https://github.com/ENCODE-DCC/wgbs-pipeline/blob/dev/wgbs_pipeline/calculate_average_coverage.py,
whose JSON layout this keeps.
"""

import argparse
import json
import subprocess


def numeric(value):
    """Convert a samtools stats value to int or float where it parses as one.

    Parameters
    ----------
    value : str
        Value field of an ``SN`` line.

    Returns
    -------
    int or float or str
        The value as a number, or unchanged.
    """
    for kind in (int, float):
        try:
            return kind(value)
        except ValueError:
            pass
    return value


def samtools_stats(bam, threads):
    """Run samtools stats on a BAM and parse its summary numbers.

    Parameters
    ----------
    bam : str
        BAM file.
    threads : int
        Threads for samtools stats.

    Returns
    -------
    dict of str to int or float or str
        One entry per ``SN`` line, keyed by its name without the colon.
    """
    out = subprocess.run(
        ["samtools", "stats", "--threads", str(threads), bam],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    fields = [line.split("\t") for line in out.splitlines() if line.startswith("SN\t")]
    return {f[1].strip(":"): numeric(f[2]) for f in fields}


def genome_size(chromsizes):
    """Sum the sizes in a chromosome sizes file.

    Parameters
    ----------
    chromsizes : str
        Two-column chromosome sizes file.

    Returns
    -------
    int
        Genome size in bases.
    """
    with open(chromsizes) as f:
        return sum(int(line.split()[1]) for line in f)


def build_parser():
    """Build the make_average_coverage.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the make_average_coverage.py flags.
    """
    parser = argparse.ArgumentParser(description="Write the mapping QC JSON.")
    parser.add_argument("--bamfile", required=True, help="Final BAM.")
    parser.add_argument("--chromsizes", required=True, help="Chromosome sizes file.")
    parser.add_argument(
        "--threads", type=int, default=1, help="Threads for samtools stats."
    )
    parser.add_argument(
        "--gem_mapper_jsons", nargs="+", required=True, help="gem-mapper JSON reports."
    )
    parser.add_argument("--output", required=True, help="QC JSON to write.")
    return parser


def main(argv=None):
    """Write the gem-mapper reports, samtools stats and average coverage as one JSON.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    report = {}
    for path in args.gem_mapper_jsons:
        with open(path) as f:
            report[path] = json.load(f)
    stats = samtools_stats(args.bamfile, args.threads)
    coverage = (
        stats["reads mapped"] * stats["average length"] / genome_size(args.chromsizes)
    )
    # qc_utils sorted the metrics and their fields by name.
    report["average_coverage"] = {"average_coverage": coverage}
    report["samtools_stats"] = dict(sorted(stats.items()))
    with open(args.output, "w") as f:
        json.dump(report, f, indent=4)


if __name__ == "__main__":
    main()
