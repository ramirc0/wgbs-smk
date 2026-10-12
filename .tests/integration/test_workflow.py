"""Integration tests on the toy fixture."""

import gzip
import hashlib
import json
import subprocess
import sys

import pyBigWig
import pytest
from conftest import FIXTURE, PREFIX, ROOT, run_ok, snakemake, tool, workdir

sys.path.insert(0, str(ROOT / "workflow/scripts"))

from make_average_coverage import numeric

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


def test_untracked_files(run):
    """Check that the workflow declares every file it writes under results/."""
    done = snakemake(run, "--list-untracked")
    assert done.returncode == 0, done.stderr[-5000:]
    # Snakemake lists the files inside the pools/ directory output as untracked.
    untracked = [
        p
        for p in (done.stdout + done.stderr).split()
        if p.startswith("results/") and not p.startswith("results/refs/pools/")
    ]
    assert untracked == []


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
    """Check that fastq_split 1 maps the input pair straight to the final BAM.

    A copy of a part would double the BAM's disk use and write time, and a hard
    link fails Snakemake's output mtime check on SLURM.
    """
    path = workdir(tmp_path, [("fastq_split: 2 ", "fastq_split: 1 ")])
    run_ok(path, "--notemp", f"{PREFIX}/S001/S001.bam")
    assert not (path / PREFIX / "S001/split").exists()
    assert not list((path / PREFIX / "S001").glob("S001.0*.bam"))
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


def test_old_sheet_headers_name_the_new_ones(tmp_path):
    """Check that a sheet with the pre-template headers fails with a rename hint."""
    path = workdir(tmp_path)
    sheet = path / "config/samples.tsv"
    lines = sheet.read_text().splitlines(keepends=True)
    sheet.write_text("Barcode\tDataset\tFile1\tFile2\n" + "".join(lines[1:]))
    done = snakemake(path, "-n")
    assert done.returncode != 0
    assert "Barcode, Dataset, File1, File2" in done.stdout + done.stderr


@pytest.mark.parametrize(
    "name", ["cpg.txt.gz", "non_cpg.txt.gz", "cpg.bed.gz", "chg.bed.gz", "chh.bed.gz"]
)
def test_mextr_index_finds_each_site_at_its_position(run, name):
    """Check that a 1-based region query on a mextr file returns the site there."""
    path = run / PREFIX / f"S001/S001_{name}"
    with gzip.open(path, "rt") as f:
        f.readline()
        contig, start, end = f.readline().split("\t")[:3]
    hits = tool("tabix", path, f"{contig}:{end}-{end}").splitlines()
    assert [h.split("\t")[:3] for h in hits] == [[contig, start, end]]


def test_mextr_keeps_the_last_base_of_a_contig(run, tmp_path):
    """Check that mextr reports a site on the contig's last base.

    The toy reads stop short of the contig ends, so the test shortens chr1 in
    the BCF header to end at its last reported site.
    """
    sample = run / PREFIX / "S001/S001"
    last = tool("tabix", f"{sample}_non_cpg.txt.gz", "chr1").splitlines()[-1]
    end = last.split("\t")[2]
    hdr = tmp_path / "header.txt"
    hdr.write_text(
        tool("bcftools", "view", "-h", f"{sample}.bcf").replace(
            "ID=chr1,length=2000", f"ID=chr1,length={end}"
        )
    )
    bcf = tmp_path / "short.bcf"
    tool("bcftools", "reheader", "-h", hdr, "-o", bcf, f"{sample}.bcf")
    tool("bcftools", "index", bcf)
    noncpg = tmp_path / "non_cpg.txt.gz"
    tool("mextr", "--compress", "--noncpgfile", noncpg, "--threshold", "10", bcf)
    with gzip.open(noncpg, "rt") as f:
        chr1 = [line for line in f if line.startswith("chr1\t")]
    assert chr1[-1] == last + "\n"


def test_coverage_matches_one_genomecov_pass(run):
    """Check that each strand's coverage bigWig holds one bedtools genomecov pass."""
    bam = run / PREFIX / "S001/S001.bam"
    for strand, name in (("-", "neg"), ("+", "pos")):
        lines = tool("bedtools", "genomecov", "-ibam", bam, "-bg", "-strand", strand)
        expected = [
            (c, int(s), int(e), float(v))
            for c, s, e, v in (line.split("\t") for line in lines.splitlines())
        ]
        bw = pyBigWig.open(str(run / PREFIX / f"S001/S001_coverage_{name}.bw"))
        ours = [(c, *i) for c in bw.chroms() for i in bw.intervals(c) or ()]
        assert sorted(ours) == sorted(expected)


def test_qc_json_matches_one_samtools_stats_pass(run):
    """Check that the QC JSON's samtools stats equal one pass over the whole BAM."""
    sample = run / PREFIX / "S001/S001"
    lines = tool("samtools", "stats", f"{sample}.bam").splitlines()
    fields = [line.split("\t") for line in lines if line.startswith("SN\t")]
    expected = {f[1].strip(":"): numeric(f[2]) for f in fields}
    report = json.loads(sample.with_suffix(".json").read_text())
    assert report["samtools_stats"] == expected


def test_methylc_holds_every_site_once_in_order(run):
    """Check that the methylC file merges the three bedMethyl files, sorted."""
    sample = run / PREFIX / "S001/S001"
    sites = []
    for context, name in (("CHH", "chh"), ("CHG", "chg"), ("CG", "cpg")):
        with gzip.open(f"{sample}_{name}.bed.gz", "rt") as f:
            next(f)
            for line in f:
                c, s, e, _, _, strand, _, _, _, cov, *_ = line.split("\t")
                sites.append((c, int(s), int(e), context, strand, cov))
    with gzip.open(f"{sample}_methylc.bed.gz", "rt") as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    ours = [(c, int(s), int(e), t, strand, cov) for c, s, e, t, _, strand, cov in rows]
    assert ours == sorted(sites, key=lambda site: site[:2])


@pytest.mark.parametrize("name", ["S001", "S002"])
def test_every_sample_finishes(run, name):
    """Check that each sample gets its methylC file."""
    assert (run / PREFIX / name / f"{name}_methylc.bed.gz").stat().st_size > 0
