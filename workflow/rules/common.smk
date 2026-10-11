from pathlib import Path

import polars as pl
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

SPLIT = [f"{part:03d}" for part in range(1, config["fastq_split"] + 1)]


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: sample_id, dataset, fastq_1, fastq_2; see config/samples.example.tsv."
    )

# All-string read so empty cells are "" (absent), not a type error.
_manifest = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")


def _column(name):
    """Map each sample to its value in sheet column `name`.

    Parameters
    ----------
    name : str
        Sample sheet column.

    Returns
    -------
    dict of str to str or None
        Value per `sample_id`, with empty cells as None. Empty when the sheet
        has no such column.
    """
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


SAMPLES = _manifest.get_column("sample_id").to_list()
DATASET_OF = _column("dataset")
FASTQ_1_OF = _column("fastq_1")
FASTQ_2_OF = _column("fastq_2")

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
    """Return the gem-mapper `@RG` header line of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str
        Tab-separated line with `ID` and `PU` set to the sample's dataset,
        `BC` to the sample ID, and an empty `SM`. The tabs are real characters,
        not `\\t` escapes.
    """
    dataset = DATASET_OF[sample].strip()
    return f"@RG\tID:{dataset}\tSM:\tBC:{sample}\tPU:{dataset}"


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
