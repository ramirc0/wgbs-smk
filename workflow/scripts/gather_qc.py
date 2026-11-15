"""Collect every sample's QC into one table, a row per sample, for the run's report."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl
from qc_outputs import (
    conversion_pct,
    depth_ratio,
    mean_methylation_pct,
    read_call_reports,
    read_mapping_qc,
    sites_at_depth,
)


def build_parser():
    """Build the gather_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the gather_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Collect each sample's QC metrics.")
    p.add_argument("--run-dir", required=True, help="Run directory of sample dirs.")
    p.add_argument("--samples", required=True, nargs="+", help="Sample IDs.")
    p.add_argument("--pools", required=True, help="Contig pool dir with pools.txt.")
    p.add_argument("--chromsizes", required=True, help="Chromosome sizes file.")
    p.add_argument("--threads", type=int, default=1, help="Samples read at once.")
    p.add_argument("--output", required=True, help="TSV to write.")
    return p


def sample_metrics(run, sample, pools, sizes):
    """Collect one sample's metrics, keyed by the ids in ``workflow/report/qc.yaml``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    sample : str
        Sample ID.
    pools : list of str
        Contig pool names, one bs_call report each.
    sizes : dict of str to int
        Length per contig.

    Returns
    -------
    dict
        One row of the metrics table.
    """
    d = run / sample
    mapping = read_mapping_qc(d / f"align/{sample}_qc.json")
    calls = read_call_reports([d / f"call/{sample}_{pool}.json" for pool in pools])
    stats = mapping["stats"]
    return {
        "sample_id": sample,
        "read_pairs": stats["raw total sequences"] // 2,
        "read_length": stats["maximum length"],
        "mapped_pct": 100 * stats["reads mapped"] / stats["raw total sequences"],
        "properly_paired_pct": stats["percentage of properly paired reads (%)"],
        "mq0_pct": 100 * stats["reads MQ0"] / stats["reads mapped"],
        "error_rate_pct": 100 * stats["error rate"],
        "insert_size": stats["insert size average"],
        "average_coverage": mapping["coverage"],
        "cpgs_10x": sites_at_depth(calls["depths"], 10),
        "passed_cpgs": calls["passed_cpgs"],
        "conversion_pct": conversion_pct(calls["profile"]),
        "cpg_methylation_pct": mean_methylation_pct(calls["methylation"]),
        "chrx_ratio": depth_ratio(mapping["contig_reads"], sizes, "X"),
        "chry_ratio": depth_ratio(mapping["contig_reads"], sizes, "Y"),
    }


def main(argv=None):
    """Write the run's metrics table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    run = Path(args.run_dir)
    pools = Path(args.pools, "pools.txt").read_text().split()
    with open(args.chromsizes) as f:
        sizes = {name: int(size) for name, size, *_ in map(str.split, f)}
    with ThreadPoolExecutor(args.threads) as pool:
        rows = list(
            pool.map(lambda s: sample_metrics(run, s, pools, sizes), args.samples)
        )
    pl.DataFrame(rows).write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
