"""Compare the integration run's CpG calls with ENCODE4's for the same replicates.

ENCODE4 v1.1.6 runs gemBS on GRCh38 with its own reference and settings, so the
calls are close, not equal. Run from the repository root after the integration
run:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/compare_encode.py
"""

import argparse
from pathlib import Path

import polars as pl

ENCODE_DATA = Path("/data/projects/encode/data")

# sample_id -> ENCODE4 v1.1.6 GRCh38 "methylation state at CpG" file.
ENCODE_CPG = {
    "ENCSR765JPCrep1": "ENCSR765JPC/ENCFF660IHA.bed.gz",
    "ENCSR765JPCrep2": "ENCSR765JPC/ENCFF328NMN.bed.gz",
    "ENCSR890UQOrep1": "ENCSR890UQO/ENCFF570TIL.bed.gz",
    "ENCSR890UQOrep2": "ENCSR890UQO/ENCFF614QHA.bed.gz",
}


def build_parser():
    """Build the compare_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the compare_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Compare CpG calls with ENCODE's.")
    p.add_argument(
        "--run",
        default="results/integration/encode-integration",
        help="Run directory.",
    )
    p.add_argument(
        "--min-coverage",
        type=int,
        default=10,
        help="Coverage both calls need for the methylation comparison.",
    )
    return p


def read_cpg(path):
    """Read a bedMethyl CpG file.

    Parameters
    ----------
    path : pathlib.Path
        Gzipped bedMethyl file, with or without a `track` line.

    Returns
    -------
    polars.DataFrame
        Columns `chrom`, `start`, `strand`, `coverage` and `percent`.
    """
    return pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        comment_prefix="track",
        columns=[0, 1, 5, 9, 10],
        new_columns=["chrom", "start", "strand", "coverage", "percent"],
        schema_overrides={"coverage": pl.Int32, "percent": pl.Float32},
    )


def compare(ours, encode, min_coverage):
    """Summarize the agreement of two CpG call sets.

    Parameters
    ----------
    ours, encode : polars.DataFrame
        Calls from `read_cpg`.
    min_coverage : int
        Coverage both calls need for the methylation statistics.

    Returns
    -------
    dict of str to float
        Site counts, the shared fraction of each set, coverage and methylation
        correlations, the mean absolute methylation difference, and each set's
        coverage-weighted mean methylation.
    """
    key = ["chrom", "start", "strand"]
    both = ours.join(encode, on=key, suffix="_encode")
    deep = both.filter(
        (pl.col("coverage") >= min_coverage)
        & (pl.col("coverage_encode") >= min_coverage)
    )
    stats = deep.select(
        pl.corr("percent", "percent_encode").alias("r_meth"),
        (pl.col("percent") - pl.col("percent_encode")).abs().mean().alias("mad_meth"),
    ).row(0, named=True)

    def weighted(df):
        cov = pl.col("coverage")
        return df.select((cov * pl.col("percent")).sum() / cov.sum()).item()

    return {
        "sites_ours": ours.height,
        "sites_encode": encode.height,
        "shared_of_ours": both.height / ours.height,
        "shared_of_encode": both.height / encode.height,
        "r_coverage": both.select(pl.corr("coverage", "coverage_encode")).item(),
        f"sites_cov{min_coverage}": deep.height,
        "r_meth": stats["r_meth"],
        "mad_meth": stats["mad_meth"],
        "mean_meth_ours": weighted(ours),
        "mean_meth_encode": weighted(encode),
    }


def main():
    """Print one row of CpG agreement statistics per finished replicate."""
    args = build_parser().parse_args()
    rows = []
    for sample, encode_file in ENCODE_CPG.items():
        ours = Path(args.run) / sample / f"{sample}_cpg.bed.gz"
        if not ours.exists():
            print(f"{sample}: no {ours}, skipped")
            continue
        stats = compare(
            read_cpg(ours), read_cpg(ENCODE_DATA / encode_file), args.min_coverage
        )
        rows.append({"sample": sample, **stats})
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200, float_precision=4):
        print(pl.DataFrame(rows))


if __name__ == "__main__":
    main()
