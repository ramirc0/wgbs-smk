"""Write each sample's report figures, from the QC files the pipeline keeps."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import plotly.graph_objects as go
import polars as pl
from qc_outputs import read_call_reports, read_mapping_qc
from qc_report import insert_size_figure, write_sample_panels

# Deepest CpG depth drawn. Deeper sites are rare and stretch the axis.
MAX_DEPTH = 100


def build_parser():
    """Build the qc_figures.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc_figures.py flags.
    """
    p = argparse.ArgumentParser(description="Write each sample's report figures.")
    p.add_argument("--run-dir", required=True, help="Run directory of sample dirs.")
    p.add_argument("--report-dir", required=True, help="Directory of the report.")
    p.add_argument("--samples", required=True, nargs="+", help="Sample IDs.")
    p.add_argument("--pools", required=True, help="Contig pool dir with pools.txt.")
    p.add_argument("--threads", type=int, default=1, help="Samples drawn at once.")
    return p


def non_conversion_figure(profile):
    """Draw non-conversion at non-CpG sites by position in the read.

    Parameters
    ----------
    profile : list of list of int
        ``profile`` of `qc_outputs.read_call_reports`.

    Returns
    -------
    plotly.graph_objects.Figure
        ``100 * a / (a + b)`` per read position, as gemBS plots it. Rises at
        the read ends point to end repair or adapter bases that need trimming.
    """
    points = [
        (position, 100 * kept / (kept + converted))
        for position, (kept, converted, *_) in enumerate(profile, start=1)
        if kept + converted
    ]
    fig = go.Figure(
        go.Scatter(
            x=[p for p, _ in points],
            y=[v for _, v in points],
            mode="lines",
            hovertemplate="Position %{x}<br>%{y:.2f}% not converted<extra></extra>",
        )
    )
    fig.update_layout(
        xaxis_title="Position in read",
        yaxis_title="Non-conversion at non-CpG sites (%)",
        yaxis_rangemode="tozero",
        showlegend=False,
    )
    return fig


def methylation_figure(histogram):
    """Draw the methylation of the CpGs that pass bs_call's filters.

    Parameters
    ----------
    histogram : list of float
        ``methylation`` of `qc_outputs.read_call_reports`, one bin per percent.

    Returns
    -------
    plotly.graph_objects.Figure
        Share of CpGs per methylation percent.
    """
    total = sum(histogram) or 1
    fig = go.Figure(
        go.Bar(
            x=list(range(len(histogram))),
            y=[100 * count / total for count in histogram],
            hovertemplate="%{x}% methylated<br>%{y:.2f}% of CpGs<extra></extra>",
        )
    )
    fig.update_layout(
        xaxis_title="CpG methylation (%)", yaxis_title="CpGs (%)", showlegend=False
    )
    return fig


def depth_figure(depths):
    """Draw the read depth of reference CpGs.

    Parameters
    ----------
    depths : dict of int to int
        ``depths`` of `qc_outputs.read_call_reports`.

    Returns
    -------
    plotly.graph_objects.Figure
        CpGs per depth up to `MAX_DEPTH`, with 10X marked.
    """
    shown = {d: n for d, n in depths.items() if d <= MAX_DEPTH}
    fig = go.Figure(
        go.Bar(
            x=list(shown),
            y=list(shown.values()),
            hovertemplate="%{x}X<br>%{y:,} CpG positions<extra></extra>",
        )
    )
    fig.add_vline(x=10, line_dash="dot", annotation_text="10X")
    fig.update_layout(
        xaxis_title="Read depth", yaxis_title="CpG positions", showlegend=False
    )
    return fig


def write_panels(sample, run, report_dir, pools):
    """Write one sample's ``samples/<sample>.js``.

    Parameters
    ----------
    sample : str
        Sample ID.
    run : pathlib.Path
        Run directory.
    report_dir : pathlib.Path
        Directory of the report HTML.
    pools : list of str
        Contig pool names, one bs_call report each.
    """
    d = run / sample
    mapping = read_mapping_qc(d / f"align/{sample}_qc.json")
    calls = read_call_reports([d / f"call/{sample}_{pool}.json" for pool in pools])
    panels = []
    if mapping["template_lengths"]:
        lengths = pl.DataFrame(
            {
                "insert_size": list(mapping["template_lengths"]),
                "pairs": list(mapping["template_lengths"].values()),
            }
        )
        panels.append(("Template lengths", insert_size_figure(lengths)))
    panels += [
        ("Non-conversion by read position", non_conversion_figure(calls["profile"])),
        ("CpG methylation", methylation_figure(calls["methylation"])),
        ("CpG read depth", depth_figure(calls["depths"])),
    ]
    write_sample_panels(report_dir / "samples" / f"{sample}.js", sample, panels)


def main(argv=None):
    """Write every sample's figures.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    draw = partial(
        write_panels,
        run=Path(args.run_dir),
        report_dir=Path(args.report_dir),
        pools=Path(args.pools, "pools.txt").read_text().split(),
    )
    with ProcessPoolExecutor(args.threads) as pool:
        list(pool.map(draw, args.samples))


if __name__ == "__main__":
    main()
