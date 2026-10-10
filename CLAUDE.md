# CLAUDE.md

WGBS Snakemake workflow on the snakemake-template conventions. README.md holds
the layout, idioms, conventions, and known issues. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run test                                # full run in .tests/integration
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` is the rule env. Rule dependencies go there, never in the root.
- gem3-mapper and gemBS-rs are `pixi-build` source packages in `workflow/envs/recipes/`. They compile on first install.
- Rules use `locked=True`. After a manifest or recipe edit, run `pixi lock` in `workflow/envs/`.
- `SOFTWARE_ENV` resolves `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Gotchas

- Shell blocks MUST call scripts as `{workflow.basedir:q}/scripts/x.py`. A cwd-relative path breaks `--directory`.
- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- The read group needs real tab characters. Build it in `read_group()`, never inside an `r"""` shell block.
- gembs_extract copies back everything in its staging dir. Its `tmpdir` subdir MUST stay apart from the `calls` one.
- `seqkit split2 -f` wipes its output dir. Keep the split FASTQs in their own `split/` dir.
- gem3-mapper is built with `-march=x86-64-v4`. It needs AVX-512 nodes.
- Scripts pin polars 1.x. They use the 1.x `read_csv` API.
- `/` on this host is nearly full. Point `TMPDIR` at a `/home` path for builds and test runs.
- Known issues in README.md are kept on purpose for parity with the `pre-template` tag. Fix one only with a test that shows it.

## Verifying a change

`pixi run test` runs every rule on the toy data in `.tests/integration/`: every contig in `CONTIGS` plus one alt contig, and two samples of directional bisulfite reads with 80% CpG methylation.
New rules or inputs MUST get matching test inputs there.

To check faithfulness, run the same inputs through a `pre-template` worktree and compare the normalized BAM, BCF, bed, and bigWig contents.
