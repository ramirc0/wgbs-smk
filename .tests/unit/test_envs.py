"""Checks on the installed pixi environments."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
ENVS = [ROOT / ".pixi/envs/default", ROOT / "workflow/envs/.pixi/envs/default"]


@pytest.mark.parametrize("env", ENVS, ids=["launcher", "rule"])
def test_no_wheel_was_built_on_this_host(env):
    """Check that every compiled PyPI package came as a portable wheel.

    A wheel built here takes the login shell's CFLAGS, which target AVX-512,
    and crashes on the broadwell nodes. Conda packages are built elsewhere.
    """
    local = []
    for wheel in env.glob("lib/python*/site-packages/*.dist-info/WHEEL"):
        if (wheel.parent / "INSTALLER").read_text().strip() == "conda":
            continue
        tags = [
            line.split()[1]
            for line in wheel.read_text().splitlines()
            if line.startswith("Tag:")
        ]
        if any(t.endswith("linux_x86_64") and "manylinux" not in t for t in tags):
            local.append(wheel.parent.name)
    assert local == []
