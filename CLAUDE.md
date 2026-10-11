# CLAUDE.md

WGBS Snakemake workflow on the snakemake-template conventions. README.md holds
the layout, idioms, conventions, and known issues. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run test                                # pytest over .tests; toy runs in scratch copies of .tests/integration
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local
pixi run snakemake --profile profiles/slurm --configfile .tests/validation/config.yaml  # ENCODE validation run
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` is the rule env. Rule dependencies go there, never in the root.
- gem3-mapper and gemBS-rs are `pixi-build` source packages in `workflow/envs/recipes/`. They compile on first install.
- The build backend MUST stay pinned. The lock does not record it, so a new release makes every job rebuild both packages from GitHub. Compute nodes may not resolve github.com.
- A checkout or edit in `workflow/envs/` also triggers a rebuild. Run `pixi install --locked` in `workflow/envs/` before launching jobs.
- Rules use `locked=True`. After a manifest or recipe edit, run `pixi lock` in `workflow/envs/`.
- `SOFTWARE_ENV` resolves `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Gotchas

- Shell blocks MUST call scripts as `{workflow.basedir:q}/scripts/x.py`. A cwd-relative path breaks `--directory`.
- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- The read group needs real tab characters. Build it in `read_group()`, never inside an `r"""` shell block.
- mextr reads its outputs back while writing them. On CephFS the reads come back short, so `gembs_extract` stages them in `tmpdir`. The toy tests run on local disk and cannot catch this.
- `seqkit split2 -f` wipes its output dir. Keep the split FASTQs in their own `split/` dir.
- gem3-mapper is built with `-march=x86-64-v4`. It needs AVX-512 nodes.
- The login shell exports `CFLAGS=-O3 -march=x86-64-v4`. A PyPI package with no wheel for our python gets built with it and crashes on broadwell nodes. Take such packages from conda-forge. `.tests/unit/test_envs.py` checks this.
- The rule env pins polars 2.x. The scripts read mextr's headerless bedMethyl by column index.
- Point `TMPDIR` at a roomy local path, such as `/tmp/$USER`, for builds and test runs. Do not set it for a SLURM launch: jobs inherit it as `resources.tmpdir`.
- gem-mapper learns insert sizes per thread, so alignments can shift with the thread count. Compare runs with different threads or `fastq_split` at the methylation level.
- Fix a bug only after a test that shows it. Known issues in README.md come from gemBS-rs.

## Verifying a change

`pixi run test` runs every rule on the toy data in `.tests/integration/`: 2 kb of chr1-22, X, Y and one alt contig, a 1 kb chrM without reads, and two samples of directional bisulfite reads with 80% CpG methylation.
New rules or inputs MUST get matching test inputs there.
For a big change, run the validation on SLURM, then `.tests/validation/compare_encode.py`.

To check faithfulness, snapshot the toy outputs before a change, or run a `post-template` worktree, and compare the normalized BAM, BCF, bed, and bigWig contents.
