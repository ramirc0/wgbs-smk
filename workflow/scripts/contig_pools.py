"""Group contigs into bs_call pools the way gemBS does, one BED per pool.

A contig of at least `limit` bases is its own pool, named after it. The smaller
contigs are spread over ceil(total / limit) pools named `Pool@1`, `Pool@2`, ...:
longest first, each into the pool with the fewest bases so far
(gemBS-rs `config/contig.rs`). Each BED lists its contigs in reference order,
and `pools.txt` lists the pools in the order of their first contig, so the
concatenated BCF follows the reference.
"""

import argparse
import heapq
from pathlib import Path


def read_sizes(path):
    """Read a chromosome sizes file.

    Parameters
    ----------
    path : str
        Two-column file of contig names and lengths, in reference order.

    Returns
    -------
    list of tuple of (str, int)
        Name and length per contig.
    """
    with open(path) as f:
        return [(name, int(size)) for name, size in (line.split()[:2] for line in f)]


def contig_pools(sizes, limit):
    """Group contigs into pools.

    Parameters
    ----------
    sizes : list of tuple of (str, int)
        Name and length per contig, in reference order.
    limit : int
        Contigs at least this long get their own pool.

    Returns
    -------
    dict of str to list of int
        Indices into `sizes` per pool, each list in reference order. Pools are
        ordered by their first contig.
    """
    pools = {name: [i] for i, (name, size) in enumerate(sizes) if size >= limit}
    small = [i for i, (_, size) in enumerate(sizes) if size < limit]
    total = sum(sizes[i][1] for i in small)
    heap = [(0, n) for n in range(-(-total // limit))]
    members = [[] for _ in heap]
    for i in sorted(small, key=lambda i: -sizes[i][1]):
        bases, n = heapq.heappop(heap)
        members[n].append(i)
        heapq.heappush(heap, (bases + sizes[i][1], n))
    pools |= {f"Pool@{n + 1}": sorted(m) for n, m in enumerate(members)}
    return dict(sorted(pools.items(), key=lambda item: item[1][0]))


def build_parser():
    """Build the contig_pools.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the contig_pools.py flags.
    """
    parser = argparse.ArgumentParser(
        description="Write one contig BED per bs_call pool."
    )
    parser.add_argument("--chromsizes", required=True, help="Chromosome sizes file.")
    parser.add_argument(
        "--limit", type=int, required=True, help="Contigs this long get their own pool."
    )
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="Output directory."
    )
    return parser


def main(argv=None):
    """Write `<pool>.bed` per pool and `pools.txt` into the output directory.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    sizes = read_sizes(args.chromsizes)
    pools = contig_pools(sizes, args.limit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for pool, indices in pools.items():
        rows = (f"{sizes[i][0]}\t0\t{sizes[i][1]}\n" for i in indices)
        (args.output_dir / f"{pool}.bed").write_text("".join(rows))
    (args.output_dir / "pools.txt").write_text("".join(f"{p}\n" for p in pools))


if __name__ == "__main__":
    main()
