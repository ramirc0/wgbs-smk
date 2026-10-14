import csv
from pathlib import Path

from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run
REFDIR = f"{RESULTS}/refs"              # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
REFERENCE = config["references"]["reference"]

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


# The pixi plugin installs the env inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install it once before any job starts.
onstart:
    shell("pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked")

SPLIT = [f"{part:03d}" for part in range(1, config["fastq_split"] + 1)]

ALIGNER = config.get("aligner", "gem3")
if ALIGNER not in ("gem3", "bwa-mem3"):
    raise WorkflowError(f"aligner is {ALIGNER!r}. It must be gem3 or bwa-mem3.")


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: sample_id, dataset, fastq_1, fastq_2; see config/samples.example.tsv."
    )

with _SHEET.open(newline="") as f:
    _rows = csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")
    _missing = {"sample_id", "dataset", "fastq_1", "fastq_2"} - set(_rows.fieldnames or ())
    if _missing:
        raise WorkflowError(
            f"Sample sheet {_SHEET} lacks the columns {sorted(_missing)}.\n"
            "Sheets from before the template port name them Barcode, Dataset, File1, File2."
        )
    SHEET = {row["sample_id"]: row for row in _rows}

SAMPLES = list(SHEET)

if not SAMPLES:
    raise WorkflowError(f"Sample sheet {_SHEET} has no rows.")


wildcard_constraints:
    sample=r"[A-Za-z0-9]+",
    part=r"\d{3}",
    pool=r"[^/]+",


def prefix(sample):
    """Return the output prefix of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<sample>/<sample>`.
    """
    return f"{OUTDIR}/{sample}/{sample}"


def read_group(sample):
    """Return the aligner's `@RG` header line of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str
        Line with `ID` and `PU` set to the sample's dataset, `BC` to the sample
        ID, and an empty `SM`. For gem-mapper, real tab characters separate the
        fields. bwa-mem3 rejects those, so its line has `\\t` escapes.
    """
    dataset = SHEET[sample]["dataset"].strip()
    line = f"@RG\tID:{dataset}\tSM:\tBC:{sample}\tPU:{dataset}"
    return line if ALIGNER == "gem3" else line.replace("\t", "\\t")


def dict_flags():
    """Return the samtools dict flags for the optional assembly and species.

    Returns
    -------
    list of str
        Flag tokens. Interpolate them with `:q`.
    """
    flags = []
    if config.get("assembly", ""):
        flags += ["-a", str(config["assembly"])]
    if config.get("species", ""):
        flags += ["-s", str(config["species"])]
    return flags


def call_flags():
    """Return the bs_call flags from `config["calling"]`, in gemBS's order.

    Returns
    -------
    list of str
        Flag tokens. Interpolate them with `:q`.
    """
    c = config["calling"]
    return [
        "--conversion", str(c["conversion"]),
        "--left-trim", str(c["left_trim"]),
        "--right-trim", str(c["right_trim"]),
        "--reference-bias", str(c["reference_bias"]),
        "--mapq-threshold", str(c["mapq_threshold"]),
        "--bq-threshold", str(c["qual_threshold"]),
    ]


def extract_flags():
    """Return the optional mextr flags from `config["extract"]`.

    Returns
    -------
    list of str
        Flag tokens. Interpolate them with `:q`.
    """
    e = config["extract"]
    flags = []
    if e["phred_threshold"]:
        flags += ["--threshold", str(e["phred_threshold"])]
    if e["strand_specific"]:
        flags += ["--mode", "strand-specific"]
    if e["bw_strand_specific"]:
        flags += ["--bw-mode", "strand-specific"]
    return flags
