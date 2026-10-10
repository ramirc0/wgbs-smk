"""Add samtools stats and average coverage to the gem-mapper QC JSON.

Coverage = (aligned read count * read length) / genome size. The genome size is
the sum of the chromosome sizes. The read length is the samtools stats average
length. Adapted from
https://github.com/ENCODE-DCC/wgbs-pipeline/blob/dev/wgbs_pipeline/calculate_average_coverage.py.
"""

import argparse
import json
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Any

import pysam
from qc_utils import QCMetric, QCMetricRecord
from qc_utils.parsers import parse_samtools_stats


def get_samtools_stats(bam_path: str, threads: int) -> dict[str, Any]:
    """Run samtools stats on a BAM and parse the summary numbers.

    The temporary file is flushed before the parser reads it. See
    https://stackoverflow.com/questions/46004774/python-namedtemporaryfile-appears-empty-even-after-data-is-written.

    Parameters
    ----------
    bam_path : str
        BAM file.
    threads : int
        Threads for samtools stats.

    Returns
    -------
    dict of str to Any
        Parsed samtools stats summary.
    """
    samtools_stats_stdout = pysam.stats("--threads", str(threads), bam_path)
    with tempfile.NamedTemporaryFile(mode="w") as samtools_stats_file:
        samtools_stats_file.write(samtools_stats_stdout)
        samtools_stats_file.flush()
        samtools_stats = parse_samtools_stats(samtools_stats_file.name)
    return samtools_stats


def calculate_genome_size(chrom_sizes_path: str) -> int:
    """Sum the sizes in a chromosome sizes file.

    Parameters
    ----------
    chrom_sizes_path : str
        Two-column chromosome sizes file.

    Returns
    -------
    int
        Genome size in bases.
    """
    size = 0
    with open(chrom_sizes_path) as f:
        for line in f:
            size += int(line.split()[1])
    return size


def calculate_average_coverage(
    genome_size: int, aligned_read_count: int, read_length: int
) -> dict[str, float]:
    """Compute average coverage from aligned reads and genome size.

    Parameters
    ----------
    genome_size : int
        Genome size in bases.
    aligned_read_count : int
        Mapped reads.
    read_length : int
        Average read length.

    Returns
    -------
    dict of str to float
        `{"average_coverage": coverage}`.
    """
    average_coverage = (aligned_read_count * read_length) / genome_size
    return {"average_coverage": average_coverage}


def make_qc_record(qcs: list[tuple[str, dict[str, Any]]]) -> QCMetricRecord:
    """Collect named metric dicts into a QC record.

    Parameters
    ----------
    qcs : list of tuple of (str, dict)
        Metric name and values.

    Returns
    -------
    qc_utils.QCMetricRecord
        Record holding one metric per entry.
    """
    qc_record = QCMetricRecord()
    for name, qc in qcs:
        metric = QCMetric(name, qc)
        qc_record.add(metric)
    return qc_record


def build_parser():
    """Build the make_average_coverage.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the make_average_coverage.py flags.
    """
    parser = argparse.ArgumentParser(
        description="Gather additional metrics from the alignment file"
    )
    parser.add_argument("--bamfile", help="Path to bam file")
    parser.add_argument("--chromsizes", help="Path to chromsizes file")
    parser.add_argument("--threads", help="Number of threads (samtools flagstat)")
    parser.add_argument(
        "--gem_mapper_jsons",
        nargs="+",
        help="List of input JSON files from gem3-mapper",
    )
    return parser


def main(argv=None):
    """Write `<bamfile stem>.json` with the gem-mapper reports and new metrics.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    samtools_stats = get_samtools_stats(args.bamfile, args.threads)
    genome_size = calculate_genome_size(args.chromsizes)
    average_coverage = calculate_average_coverage(
        genome_size=genome_size,
        aligned_read_count=samtools_stats["reads mapped"],
        read_length=samtools_stats["average length"],
    )
    qc_record = make_qc_record(
        [("samtools_stats", samtools_stats), ("average_coverage", average_coverage)]
    )

    existing_data = OrderedDict()
    for json_file in args.gem_mapper_jsons:
        with open(json_file) as f:
            existing_data.update(
                {json_file: json.load(f, object_pairs_hook=OrderedDict)}
            )

    existing_data.update(qc_record.to_ordered_dict())

    json_path = Path(args.bamfile).parent / Path(args.bamfile).name
    with json_path.with_suffix(".json").open("w") as f:
        json.dump(existing_data, f, indent=4)


if __name__ == "__main__":
    main()
