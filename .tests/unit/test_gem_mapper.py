"""Checks on the gem-mapper build of the rule environment."""

import gzip
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENOME = ROOT / ".tests/integration/resources/genome.fa.gz"

# gem-mapper reads each input in 32 MiB buffers. It hands a thread blocks of at
# least 4 MiB that end on a multiple of 8 lines.
BUFFER = 32 << 20
BLOCK = 4 << 20
# With 256-byte records, 8 lines are 512 bytes and every block is exactly 4 MiB.
READ = 110
NAME = "read{:026d}"


def write_pair(genome, r1, r2, split):
    """Write a FASTQ pair whose read 1 has a header across the first buffer end.

    A record of the buffer's last block gets a name `split` bytes shorter. The
    block then ends `split` bytes short of 4 MiB, at the buffer's last full
    line, and the next header's first `split` bytes fill it.

    Parameters
    ----------
    genome : str
        Sequence the reads are copied from.
    r1, r2 : pathlib.Path
        Gzipped FASTQ files to write.
    split : int
        Bytes of the header that fall before the buffer end.

    Returns
    -------
    list of str
        Read names, one per pair.
    """
    complement = str.maketrans("ACGT", "TGCA")
    record = 2 * (READ + 1) + len(NAME.format(0)) + 4
    names = [NAME.format(i) for i in range(BUFFER // record + 64)]
    short = (BUFFER - BLOCK) // record + 2
    names[short] = names[short][split:]
    with (
        gzip.open(r1, "wt", compresslevel=1) as f1,
        gzip.open(r2, "wt", compresslevel=1) as f2,
    ):
        for i, name in enumerate(names):
            start = (i * 37) % (len(genome) - 300)
            mate = genome[start + 300 - READ : start + 300].translate(complement)[::-1]
            f1.write(f"@{name}\n{genome[start : start + READ]}\n+\n{'I' * READ}\n")
            f2.write(f"@{name}\n{mate}\n+\n{'I' * READ}\n")
    return names


def test_header_across_a_read_buffer_end_is_kept(tmp_path):
    """Check that gem-mapper maps a pair whose header spans two read buffers.

    gem-mapper 882b0d8 closed a block at a buffer end when the bytes of the
    partial line there filled it to 4 MiB. It dropped those bytes and the next
    block began mid-line: "Beginning Symbol ('>' or '@') not found".
    """
    with gzip.open(GENOME, "rt") as f:
        contigs = f.read()
    (tmp_path / "genome.fa").write_text(contigs)
    genome = (
        "".join(line for line in contigs.splitlines() if not line.startswith(">"))
        .upper()
        .replace("N", "A")
    )
    names = write_pair(
        genome, tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz", split=8
    )
    subprocess.run(
        ["gem-indexer", "-i", "genome.fa", "-o", "genome", "--bisulfite-index"]
        + ["--tmp-folder", f"{tmp_path}/"],
        cwd=tmp_path,
        capture_output=True,
        check=True,
    )
    done = subprocess.run(
        ["gem-mapper", "--threads", "2", "-I", "genome.gem"]
        + ["--i1", "r1.fastq.gz", "--i2", "r2.fastq.gz", "--paired-end-alignment"]
        + ["--bisulfite-conversion", "inferred-C2T-G2A"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert done.returncode == 0, done.stderr[-2000:]
    mapped = Counter(
        line.split("\t", 1)[0]
        for line in done.stdout.splitlines()
        if not line.startswith("@")
    )
    assert mapped == Counter({name: 2 for name in names})
