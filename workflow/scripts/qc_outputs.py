"""Read the QC files this pipeline writes, for the run's report."""

import json
import re
from pathlib import Path

AUTOSOME = re.compile(r"(chr)?\d+")


def read_mapping_qc(path):
    """Read the mapping QC JSON that make_average_coverage.py writes.

    Parameters
    ----------
    path : str or pathlib.Path
        ``align/<sample>_qc.json``.

    Returns
    -------
    dict
        ``stats``, the samtools stats summary numbers; ``coverage``, the average
        coverage; ``contig_reads``, mapped reads per contig; and
        ``template_lengths``, gem-mapper's template length histogram summed over
        its reports, empty for bwa-mem3.
    """
    report = json.loads(Path(path).read_text())
    lengths = {}
    for key, value in report.items():
        if key.endswith(".map.json"):
            for length, count in value["HistTemplateLen"].items():
                lengths[int(length)] = lengths.get(int(length), 0) + count
    return {
        "stats": report["samtools_stats"],
        "coverage": report["average_coverage"]["average_coverage"],
        "contig_reads": report["contig_reads_mapped"],
        "template_lengths": dict(sorted(lengths.items())),
    }


def add_rows(total, rows):
    """Add rows of counts element-wise into a running total.

    Parameters
    ----------
    total : list
        Running total, extended where ``rows`` is longer.
    rows : list
        Numbers, or lists of numbers, one per position.

    Returns
    -------
    list
        The new total.
    """
    total = total + [None] * (len(rows) - len(total))
    for i, row in enumerate(rows):
        if total[i] is None:
            total[i] = row
        elif isinstance(row, list):
            total[i] = [a + b for a, b in zip(total[i], row, strict=True)]
        else:
            total[i] += row
    return total


def read_call_reports(paths):
    """Sum a sample's bs_call reports, one per contig pool.

    Parameters
    ----------
    paths : list of str or pathlib.Path
        ``call/<sample>_<pool>.json`` files.

    Returns
    -------
    dict
        ``profile``, the non-CpG read profile: per read position, the counts
        ``[a, b, c, d]`` of bs_call's ``NonCpGreadProfile``. ``methylation``,
        the methylation histogram of passed reference CpGs, one bin per percent.
        ``depths``, reference CpG positions per read depth. ``passed_cpgs``,
        the reference CpG positions that pass bs_call's filters. bs_call counts
        a CpG's C and G as two positions.
    """
    profile, methylation, depths, passed = [], [], {}, 0
    for path in paths:
        stats = json.loads(Path(path).read_text())["totalStats"]
        profile = add_rows(profile, stats["methylation"]["NonCpGreadProfile"])
        methylation = add_rows(methylation, stats["methylation"]["PassedRefCpg"])
        for depth, count in stats["coverage"]["RefCpG"].items():
            depths[int(depth)] = depths.get(int(depth), 0) + count
        passed += stats["RefCpG"]["Passed"]
    return {
        "profile": profile,
        "methylation": methylation,
        "depths": dict(sorted(depths.items())),
        "passed_cpgs": passed,
    }


def conversion_pct(profile):
    """Estimate the C to T conversion rate from the non-CpG read profile.

    Parameters
    ----------
    profile : list of list of int
        ``profile`` of `read_call_reports`.

    Returns
    -------
    float or None
        ``100 * b / (a + b)`` over all read positions: the complement of the
        non-conversion gemBS plots, ``a / (a + b)``. Columns ``a`` and ``b``
        hold the Cs bs_call reads as kept and converted at non-CpG sites on the
        converted strand. None without such sites.
    """
    kept = sum(row[0] for row in profile)
    converted = sum(row[1] for row in profile)
    total = kept + converted
    return 100 * converted / total if total else None


def mean_methylation_pct(histogram):
    """Average a methylation histogram whose bin i stands for i percent.

    Parameters
    ----------
    histogram : list of float
        ``methylation`` of `read_call_reports`.

    Returns
    -------
    float or None
        Mean methylation in percent, or None for an empty histogram.
    """
    total = sum(histogram)
    if not total:
        return None
    return sum(i * count for i, count in enumerate(histogram)) / total


def sites_at_depth(depths, min_depth):
    """Count the sites at or above a read depth.

    Parameters
    ----------
    depths : dict of int to int
        Sites per read depth.
    min_depth : int
        Lowest depth counted.

    Returns
    -------
    int
        Sites with at least ``min_depth`` reads.
    """
    return sum(count for depth, count in depths.items() if depth >= min_depth)


def depth_ratio(contig_reads, sizes, name):
    """Compare a sex chromosome's reads per base with the autosomes'.

    Parameters
    ----------
    contig_reads : dict of str to int
        Mapped reads per contig.
    sizes : dict of str to int
        Length per contig.
    name : str
        ``X`` or ``Y``. The contig is ``chrX`` or ``X``.

    Returns
    -------
    float or None
        About 1 for two copies per diploid genome, 0.5 for one. None when the
        reference has no such contig or no numbered autosomes.
    """
    contig = next((c for c in (f"chr{name}", name) if c in sizes), None)
    autosomes = [c for c in sizes if AUTOSOME.fullmatch(c)]
    if contig is None or not autosomes:
        return None
    autosome_depth = sum(contig_reads.get(c, 0) for c in autosomes) / sum(
        sizes[c] for c in autosomes
    )
    if not autosome_depth:
        return None
    return contig_reads.get(contig, 0) / sizes[contig] / autosome_depth
