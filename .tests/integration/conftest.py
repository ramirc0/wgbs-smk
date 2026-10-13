"""Helpers for the integration tests: scratch copies of the toy fixture and workflow runs."""

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / ".tests/integration"
PREFIX = "results/default"


def snakemake(workdir, *args):
    """Run the launcher's Snakemake on a working directory with the local profile.

    Parameters
    ----------
    workdir : pathlib.Path
        Working directory, as made by `workdir`.
    *args : str
        More Snakemake arguments, such as targets.

    Returns
    -------
    subprocess.CompletedProcess
        The finished process, with stdout and stderr captured as text.
    """
    return subprocess.run(
        ["pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"), "snakemake"]
        + ["--snakefile", str(ROOT / "workflow/Snakefile"), "--directory", str(workdir)]
        + ["--profile", "profiles/local", "--cores", "2", *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def workdir(path, config_edits=()):
    """Copy the fixture's config and resources into a scratch working directory.

    Parameters
    ----------
    path : pathlib.Path
        Directory to fill, e.g. pytest's ``tmp_path``.
    config_edits : sequence of tuple of (str, str)
        Text replacements applied to ``config/config.yaml``, each found once.

    Returns
    -------
    pathlib.Path
        The working directory.
    """
    for name in ("config", "resources"):
        shutil.copytree(FIXTURE / name, path / name)
    config = path / "config/config.yaml"
    text = config.read_text()
    for old, new in config_edits:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    config.write_text(text)
    return path


def run_ok(path, *args):
    """Run the workflow and fail the test with its log tail if it fails.

    Parameters
    ----------
    path : pathlib.Path
        Working directory.
    *args : str
        More Snakemake arguments, such as targets.
    """
    done = snakemake(path, *args)
    assert done.returncode == 0, done.stderr[-5000:]


def tool(*args):
    """Run a tool from the test environment and return its stdout.

    Parameters
    ----------
    *args : str or pathlib.Path
        Command and arguments.

    Returns
    -------
    str
        Standard output.
    """
    return subprocess.run(
        [str(a) for a in args], capture_output=True, text=True, check=True
    ).stdout


@pytest.fixture(scope="session", params=["gem3", "bwa-mem3"])
def aligner(request):
    """Return each aligner in turn.

    Returns
    -------
    str
        Value of the ``aligner`` config key.
    """
    return request.param


@pytest.fixture(scope="session")
def run(tmp_path_factory, aligner):
    """Run the whole workflow once per aligner on the toy fixture.

    Returns
    -------
    pathlib.Path
        Working directory of the run.
    """
    path = workdir(
        tmp_path_factory.mktemp(f"run-{aligner}"),
        [("aligner: gem3 ", f"aligner: {aligner} ")],
    )
    run_ok(path)
    return path
