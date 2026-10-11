"""Integration tests on the toy fixture."""

import gzip
import hashlib
import subprocess

import pytest
from conftest import FIXTURE, PREFIX, ROOT, run_ok, snakemake, tool, workdir

# bigWig files start with this magic number, little-endian.
BIGWIG_MAGIC = bytes.fromhex("26fc8f88")


def fasta_contigs(path):
    """Read a gzipped FASTA's contigs in file order.

    Parameters
    ----------
    path : pathlib.Path
        Gzipped FASTA.

    Returns
    -------
    dict of str to str
        Sequence per contig name.
    """
    contigs = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                contigs[name] = []
            else:
                contigs[name].append(line.strip())
    return {name: "".join(parts) for name, parts in contigs.items()}


def header(bam, kind):
    """Return the header lines of one kind of a BAM.

    Parameters
    ----------
    bam : pathlib.Path
        BAM file.
    kind : str
        Header record type, such as ``@RG``.

    Returns
    -------
    list of list of str
        Tab-separated fields per line.
    """
    lines = tool("samtools", "view", "-H", bam).splitlines()
    return [line.split("\t") for line in lines if line.startswith(kind + "\t")]


def test_one_read_group_per_dataset(run):
    """Check that merging the FASTQ parts keeps a single read group per dataset."""
    bam = run / PREFIX / "S001/S001.bam"
    assert [f[1] for f in header(bam, "@RG")] == ["ID:EXP1_rep1"]
    tags = tool("samtools", "view", bam).splitlines()
    assert {f for line in tags for f in line.split("\t") if f.startswith("RG:Z:")} == {
        "RG:Z:EXP1_rep1"
    }


def test_sq_lines_carry_md5_assembly_and_species(run):
    """Check that read_filter adds M5, AS and SP to the BAM's @SQ lines, as gemBS does."""
    contigs = fasta_contigs(FIXTURE / "resources/genome.fa.gz")
    sq = {f[1][3:]: f[2:] for f in header(run / PREFIX / "S001/S001.bam", "@SQ")}
    assert sq.keys() == contigs.keys()
    for name, tags in sq.items():
        md5 = hashlib.md5(contigs[name].upper().encode()).hexdigest()
        assert tags == [
            f"LN:{len(contigs[name])}",
            f"M5:{md5}",
            "AS:GRCh38",
            "SP:Homo_sapiens",
        ]


@pytest.mark.xfail(strict=True, reason="BCF and mextr md5 files name the staging path")
def test_md5_sidecars_verify(run):
    """Check that every .md5 file names its file by a path that exists."""
    sidecars = sorted((run / "results").rglob("*.md5"))
    assert sidecars
    for md5 in sidecars:
        done = subprocess.run(
            ["md5sum", "-c", md5], cwd=run, capture_output=True, text=True, check=False
        )
        assert done.returncode == 0, f"{md5.name}: {done.stdout}{done.stderr}"


def test_merged_bcf_follows_reference_order(run):
    """Check that the merged BCF lists contigs in reference order, whatever the locale."""
    body = tool("bcftools", "view", "-H", run / PREFIX / "S001/S001.bcf").splitlines()
    seen = list(dict.fromkeys(line.split("\t", 1)[0] for line in body))
    order = list(fasta_contigs(FIXTURE / "resources/genome.fa.gz"))
    assert seen == [c for c in order if c in seen]


def test_contig_beds_are_built_once_per_reference(run):
    """Check that the contig BEDs for calling live with the reference, not per sample."""
    beds = sorted((run / "results/refs").rglob("*.bed"))
    names = {bed.stem for bed in beds}
    assert len(beds) == 26
    assert "Pool@1" in names


@pytest.mark.xfail(strict=True, reason="mextr, tabix and QC files are not declared")
def test_untracked_files(run):
    """Check that the workflow declares every file it writes under results/."""
    done = snakemake(run, "--list-untracked")
    assert done.returncode == 0, done.stderr[-5000:]
    untracked = [
        p for p in (done.stdout + done.stderr).split() if p.startswith("results/")
    ]
    assert untracked == []


@pytest.mark.xfail(strict=True, reason="staging dirs stay behind")
def test_tmpdir_is_empty_after_the_run(run):
    """Check that no job leaves files or directories in tmpdir."""
    assert list((run / "scratch").rglob("*")) == []


def test_intermediates_are_removed(run):
    """Check that split FASTQs, per-part BAMs and coverage bedGraphs are temporary."""
    sample = run / PREFIX / "S001"
    leftovers = [
        *sample.glob("split/*.fastq.gz"),
        *sample.glob("S001.[0-9]*.bam"),
        *sample.glob("*.bg"),
    ]
    assert leftovers == []


def test_fastq_split_one_maps_the_pair_whole(tmp_path, run):
    """Check that fastq_split 1 maps the input pair without splitting it."""
    path = workdir(tmp_path, [("fastq_split: 2 ", "fastq_split: 1 ")])
    run_ok(path, f"{PREFIX}/S001/S001.bam")
    assert not (path / PREFIX / "S001/split").exists()
    count = [
        tool("samtools", "view", "-c", p / PREFIX / "S001/S001.bam")
        for p in (path, run)
    ]
    assert count[0] == count[1]


def test_fastq_split_into_three_parts(tmp_path):
    """Check that fastq_split other than 2 maps every part."""
    path = workdir(tmp_path, [("fastq_split: 2 ", "fastq_split: 3 ")])
    run_ok(path, f"{PREFIX}/S001/S001.bam")
    assert header(path / PREFIX / "S001/S001.bam", "@RG")


def test_reference_named_reference_fasta(tmp_path):
    """Check that a plain reference called reference.fasta is indexed."""
    path = workdir(tmp_path, [("resources/genome.fa.gz", "resources/reference.fasta")])
    with gzip.open(path / "resources/genome.fa.gz", "rb") as f:
        (path / "resources/reference.fasta").write_bytes(f.read())
    run_ok(path, "results/refs/reference.fasta.fai")


def test_contigs_without_chr_names(tmp_path):
    """Check that a reference without chr1-22, X and Y calls every contig."""
    path = workdir(tmp_path)
    fasta = path / "resources/genome.fa.gz"
    with gzip.open(fasta, "rt") as f:
        text = f.read().replace(">chr", ">ctg")
    with gzip.open(fasta, "wt") as f:
        f.write(text)
    run_ok(path)
    body = tool("bcftools", "view", "-H", path / PREFIX / "S001/S001.bcf").splitlines()
    assert {line.split("\t", 1)[0] for line in body} >= {"ctg1", "ctgX"}


@pytest.mark.xfail(strict=True, reason="extract_signal writes an empty file")
def test_empty_strand_gives_a_valid_bigwig(tmp_path):
    """Check that a context and strand without sites still gets a bigWig."""
    track = 'track name = "S" description="S" visibility=2 itemRgb="On"\n'
    site = 'chr1\t32\t33\t"S"\t1\t+\t32\t33\t255,0,0\t1\t100\tCG\tCG\t2\n'
    beds = []
    for context in ("chg", "chh", "cpg"):
        bed = tmp_path / f"S_{context}.bed.gz"
        with gzip.open(bed, "wt") as f:
            f.write(track + site)
        beds.append(bed)
    sizes = tmp_path / "chromsizes.tsv"
    sizes.write_text("chr1\t2000\n")
    subprocess.run(
        ["python", ROOT / "workflow/scripts/extract_signal.py", "--bed-gz", *beds]
        + ["--chromsizes", sizes, "--output-dir", tmp_path, "--barcode", "S"],
        check=True,
    )
    for context in ("chg", "chh", "cpg"):
        for strand in ("pos", "neg"):
            assert (tmp_path / f"S_{context}_{strand}.bw").read_bytes()[
                :4
            ] == BIGWIG_MAGIC


@pytest.mark.parametrize("name", ["S001", "S002"])
def test_every_sample_finishes(run, name):
    """Check that each sample gets its methylC file."""
    assert (run / PREFIX / name / f"{name}_methylc.bed.gz").stat().st_size > 0
