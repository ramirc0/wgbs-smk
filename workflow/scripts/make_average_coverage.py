"""Write the mapping QC JSON: the gem-mapper reports, samtools stats and average coverage.

Coverage = (mapped reads * average read length) / genome size. The genome size is
the sum of the chromosome sizes. Adapted from
https://github.com/ENCODE-DCC/wgbs-pipeline/blob/dev/wgbs_pipeline/calculate_average_coverage.py,
whose JSON layout this keeps.

samtools stats runs on one thread, about an hour on a 145 GB BAM. So it runs per
contig in parallel, and `merge_stats` rebuilds the whole-BAM summary numbers
the way samtools 1.24's `output_stats` derives them.
"""

import argparse
import json
import math
import struct
import subprocess
from concurrent.futures import ThreadPoolExecutor

# Summary numbers that add up across contigs, and those that take the maximum.
SUMMED = (
    "raw total sequences",
    "filtered sequences",
    "sequences",
    "1st fragments",
    "last fragments",
    "reads mapped",
    "reads mapped and paired",
    "reads unmapped",
    "reads properly paired",
    "reads paired",
    "reads duplicated",
    "reads MQ0",
    "reads QC failed",
    "non-primary alignments",
    "supplementary alignments",
    "total length",
    "total first fragment length",
    "total last fragment length",
    "bases mapped",
    "bases mapped (cigar)",
    "bases trimmed",
    "bases duplicated",
    "mismatches",
)
LONGEST = (
    "maximum length",
    "maximum first fragment length",
    "maximum last fragment length",
)


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


def f32(value):
    """Round a number to the nearest C float.

    Parameters
    ----------
    value : int or float
        Number to round.

    Returns
    -------
    float
        The C float value, as a Python float.
    """
    return struct.unpack("f", struct.pack("f", value))[0]


def f32_ratio(numerator, denominator):
    """Divide two counts as samtools does, ``(float)a / b``, or 0 for ``b == 0``.

    Parameters
    ----------
    numerator, denominator : int
        Counts.

    Returns
    -------
    float
        The C float quotient.
    """
    return f32(f32(numerator) / f32(denominator)) if denominator else 0


def region_output(bam, region, command):
    """Run a samtools command on the reads of one BAM region.

    Parameters
    ----------
    bam : str
        Indexed BAM file.
    region : str
        Contig name, or ``*`` for reads without a position.
    command : list of str
        Command that reads the region's reads as BAM from stdin.

    Returns
    -------
    str
        The command's standard output.
    """
    view = subprocess.Popen(
        ["samtools", "view", "-u", bam, region], stdout=subprocess.PIPE
    )
    tool = subprocess.Popen(
        command, stdin=view.stdout, stdout=subprocess.PIPE, text=True
    )
    view.stdout.close()
    out = tool.communicate()[0]
    if view.wait() != 0 or tool.returncode != 0:
        raise SystemExit(f"{command[:2]} failed on region {region}")
    return out


def region_stats(bam, region):
    """Run samtools stats and flagstat on one BAM region.

    Parameters
    ----------
    bam : str
        Indexed BAM file.
    region : str
        Contig name, or ``*`` for reads without a position.

    Returns
    -------
    tuple of (str, int)
        Output of samtools stats with the full insert size histogram, and the
        number of reads whose mate maps to another contig.
    """
    stats = region_output(bam, region, ["samtools", "stats", "-m", "1", "-"])
    flagstat = region_output(bam, region, ["samtools", "flagstat", "-O", "tsv", "-"])
    for line in flagstat.splitlines():
        passed, failed, name = line.split("\t")
        if name == "with mate mapped to a different chr":
            return stats, int(passed) + int(failed)
    raise SystemExit(f"flagstat has no different-chr line for region {region}")


def sections(output):
    """Split samtools stats output into its tab-separated rows by line type.

    Parameters
    ----------
    output : str
        Output of samtools stats.

    Returns
    -------
    dict of str to list of list of str
        Rows per line type, such as ``SN`` or ``IS``, without the type field.
    """
    rows = {}
    for line in output.splitlines():
        kind, _, rest = line.partition("\t")
        rows.setdefault(kind, []).append(rest.split("\t"))
    return rows


def insert_sizes(rows):
    """Summarise merged insert size rows as samtools stats does.

    Parameters
    ----------
    rows : list of list of str
        ``IS`` rows of several outputs: insert size, total, inward, outward and
        other pairs.

    Returns
    -------
    dict of str to str
        Insert size average and standard deviation over the main bulk, and the
        pairs per orientation.
    """
    counts = {}
    orientations = [0, 0, 0]
    for isize, _, *pairs in rows:
        counts[int(isize)] = counts.get(int(isize), 0) + sum(map(int, pairs))
        orientations = [a + int(b) for a, b in zip(orientations, pairs, strict=True)]
    nisize = sum(orientations)
    ibulk, bulk, average = 0, 0.0, 0.0
    for isize in range(max(counts, default=-1) + 1):
        num = counts.get(isize, 0)
        if num > 0:
            ibulk = isize + 1
        bulk += num
        average += isize * num
        # samtools stats' -m default, a C float.
        if nisize and bulk / nisize > f32(0.99):
            ibulk = isize + 1
            nisize = int(bulk)
            break
    average /= nisize or 1
    variance = 0.0
    for isize in range(1, ibulk):
        variance += (
            counts.get(isize, 0) * (isize - average) * (isize - average) / (nisize or 1)
        )
    inward, outward, other = orientations
    return {
        "insert size average": f"{average:.1f}",
        "insert size standard deviation": f"{math.sqrt(variance):.1f}",
        "inward oriented pairs": str(inward),
        "outward oriented pairs": str(outward),
        "pairs with other orientation": str(other),
    }


def merge_stats(outputs, different_chr):
    """Rebuild the whole-BAM summary numbers from per-region samtools stats.

    Parameters
    ----------
    outputs : list of str
        Outputs of samtools stats, run with ``-m 1``, on regions that together
        hold each read once.
    different_chr : int
        Reads whose mate maps to another contig, from samtools flagstat. Each
        region's stats halves its own count and rounds it down.

    Returns
    -------
    dict of str to str
        Value per ``SN`` name, as samtools stats prints it for the whole BAM.
    """
    parsed = [sections(output) for output in outputs]
    sn = [{row[0].rstrip(":"): row[1] for row in p["SN"]} for p in parsed]
    total = {name: sum(int(s[name]) for s in sn) for name in SUMMED}
    merged = {name: str(value) for name, value in total.items()}
    merged |= {name: str(max(int(s[name]) for s in sn)) for name in LONGEST}
    merged["is sorted"] = str(min(int(s["is sorted"]) for s in sn))
    sum_qual = sum(
        quality * int(count)
        for p in parsed
        for kind in ("FFQ", "LFQ")
        for row in p.get(kind, [])
        for quality, count in enumerate(row[1:])
    )
    length = total["total length"]
    ratio = f32_ratio(total["mismatches"], total["bases mapped (cigar)"])
    merged["error rate"] = f"{ratio:e}"
    merged["average length"] = f"{f32_ratio(length, total['sequences']):.0f}"
    for fragment, reads in (("first", "1st fragments"), ("last", "last fragments")):
        average = f32_ratio(total[f"total {fragment} fragment length"], total[reads])
        merged[f"average {fragment} fragment length"] = f"{average:.0f}"
    merged["average quality"] = f"{sum_qual / length if length else 0:.1f}"
    merged |= insert_sizes([row for p in parsed for row in p.get("IS", [])])
    merged["pairs on different chromosomes"] = str(different_chr // 2)
    paired = f32_ratio(100 * total["reads properly paired"], total["sequences"])
    merged["percentage of properly paired reads (%)"] = f"{paired:.1f}"
    return merged


def samtools_stats(bam, chromsizes, threads):
    """Run samtools stats over a BAM, one contig per job, and merge the results.

    Parameters
    ----------
    bam : str
        Indexed BAM file.
    chromsizes : str
        Chromosome sizes file naming the BAM's contigs.
    threads : int
        Contigs to process at once.

    Returns
    -------
    dict of str to int or float or str
        One entry per ``SN`` line of a whole-BAM samtools stats run.
    """
    with open(chromsizes) as f:
        sizes = [line.split()[:2] for line in f]
    # Largest first, so the last jobs are short. "*" holds the unplaced reads.
    regions = ["*"] + [name for name, size in sorted(sizes, key=lambda s: -int(s[1]))]
    with ThreadPoolExecutor(threads) as pool:
        results = list(pool.map(lambda region: region_stats(bam, region), regions))
    merged = merge_stats([r[0] for r in results], sum(r[1] for r in results))
    return {name: numeric(value) for name, value in merged.items()}


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
        "--threads",
        type=int,
        default=1,
        help="Contigs to run samtools stats on at once.",
    )
    parser.add_argument(
        "--gem_mapper_jsons",
        nargs="*",
        default=[],
        help="gem-mapper JSON reports. bwa-mem3 writes none.",
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
    stats = samtools_stats(args.bamfile, args.chromsizes, args.threads)
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
