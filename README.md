# WGBS Snakemake Pipeline

A Snakemake workflow for whole genome bisulfite sequencing (WGBS). It maps
paired-end reads with [GEM3](https://github.com/smarco/gem3-mapper), calls
methylation with [gemBS-rs](https://github.com/heathsc/gemBS-rs) `bs_call`, and
extracts bedMethyl, methylC and bigWig tracks with `mextr`. It follows the
`snakemake-template` layout and conventions.

## Layout

```
.tests/
  integration/           toy genome + reads, run by `pixi run test`
  validation/            four ENCODE replicates on GRCh38, compared with ENCODE4's calls
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per sample: sample_id, dataset, fastq_1, fastq_2
workflow/
  Snakefile              includes the rules + rule all
  rules/
    common.smk           config, sample sheet, path constants, helpers
    indexing.smk         reference FASTA, chromsizes, contig metadata, GEM index
    mapping.smk          split, map, merge, filter, strand coverage
    call.smk             contig pools, per-pool bs_call, merged BCF
    extract.smk          mextr, strand bigWigs, methylC
  scripts/               argparse scripts, each with build_parser() + main()
  envs/
    pixi.toml            rule env (pixi workspace)
    recipes/             gem3-mapper and gemBS-rs, built from source
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + test and lint tasks
pyproject.toml           ruff + numpydoc config
```

## Idioms

- **Outputs grouped by `run_id`** under `results/<run_id>/<sample>/`. Logs and
  benchmarks mirror it (`logs/<run_id>/`, `benchmarks/<run_id>/`). The
  reference artifacts are run-independent and go in `results/refs/`.
- **Sample sheet** is a TSV read once in `common.smk` into `SHEET`, one row
  dict per `sample_id`.
- **Config -> flags**: helpers such as `call_flags()`, `extract_flags()` and `dict_flags()` turn
  `config` into flag lists. `read_group()` builds the gem-mapper `@RG` line.
- **Every rule** carries `log:`, `benchmark:`, and `software:`.

## Conventions

Follow the [Nextstrain Snakemake style guide][style-guide] and
[Snakemake's best-practices guide][best-practices]:

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation quoted with `:q`. Multi-value params passed as lists.
- Config mapped into `shell` via `params:` (a `lambda`), never interpolated directly.
- `config[key]` for required keys, `config.get(key, default)` for optional ones.
- `input:`/`output:` are literal path strings (no `rules.x.output` refs). All
  paths relative.
- Logs tee'd to file and terminal (`exec &> >(tee {log:q})`). The profiles set
  `--show-failed-logs`.
- Every rule has a `benchmark:`. No `message:` attribute.
- Scripts are called as `{workflow.basedir:q}/scripts/x.py` so the workflow
  runs from any `--directory`.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and
  `numpydoc lint`. Neither tool parses `.smk` files.

[style-guide]: https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html
[best-practices]: https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/slurm           # dry run
snakemake --profile profiles/slurm                 # SLURM
snakemake --profile profiles/local                 # local
pixi run test                                      # toy end-to-end tests (pytest)
```

Run from the repository root. Local runs suit small data only.

The first job builds the rule env, which compiles gem3-mapper and gemBS-rs.
Build it once before launching many jobs:

```bash
pixi install --locked --manifest-path workflow/envs/pixi.toml
```

Tune `set-resources` in the profiles from the `benchmarks/<run_id>/` TSVs after
a first run.

### Sample sheet

One row per sample. Lines starting with `#` are ignored.

- `sample_id`: unique ID, letters and digits only.
- `dataset`: dataset or run name. It fills the read group `ID` and `PU`.
- `fastq_1`, `fastq_2`: paired, gzipped FASTQ paths.

### Configuration

- `samples`, `outdir`, `run_id`, `pixi_env`: as in the template.
- `tmpdir`: node-local scratch for `samtools sort` and `gem-indexer` temporary
  files.
- `references.reference`: reference FASTA, plain or gzipped.
- `assembly`, `species` (optional): `AS` and `SP` tags of the SAM sequence
  dictionary. `species` MUST NOT contain whitespace.
- `index.sampling_rate`: gem-indexer text sampling rate. Lower values give a
  larger, faster index.
- `fastq_split`: number of parts per FASTQ pair. With 1 (the default), one job
  maps the pair whole, as gemBS maps a dataset. With more, `seqkit split2`
  splits the pair and the parts are mapped in parallel and merged.
- `calling.*`: `bs_call` output type, MAPQ and base-quality thresholds,
  reference bias, read trimming, and conversion rates
  (`under,over`, comma-delimited). `contig_pool_limit` (default 25 Mb, as in
  gemBS) sets which contigs get their own `bs_call` job.
- `extract.*`: `mextr` strand-specific bedMethyl and bigWig modes, and the
  Phred threshold.

## Tests

`pixi run test` runs every rule on the toy data in `.tests/integration/`: 2 kb
of chr1-22, X, Y and one alt contig, a 1 kb chrM without reads, and two samples
of simulated directional bisulfite reads with 80% CpG methylation.

The validation runs four ENCODE WGBS replicates end to end on full GRCh38: two
each of K562 (ENCSR765JPC) and GM12878 (ENCSR890UQO). Its inputs live on the
Weng Lab cluster, so it runs only there.

```bash
snakemake --profile profiles/slurm --configfile .tests/validation/config.yaml
srun -p 4hours -c 8 --mem=200G -t 60 \
    pixi run --manifest-path workflow/envs/pixi.toml python .tests/validation/compare_encode.py
```

Outputs go under `results/validation/`, with their own GRCh38 index.
`compare_encode.py` joins each replicate's CpG bedMethyl with ENCODE4 v1.1.6's
on chrom, start and strand. On the 2026-10-11 run, which took 13 hours, 99% of
ENCODE's CpG sites were shared, and methylation correlated at r = 0.983-0.990
where both had coverage 10 or more. ENCODE4 mapped R1 only, so its coverage is
about half of ours.

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main, the pixi
software-deployment plugin and the SLURM executor plugins. The rule env is a
second pixi workspace, `workflow/envs/`. Its `pixi.lock` pins every package.
pixi MUST be on `PATH`.

Most tools come from conda-forge and bioconda. Two are built from source by
`pixi-build`:

- gem3-mapper at `882b0d8` on its `gembs` branch, the commit gemBS-rs pins,
  compiled with `-O3 -march=x86-64-v4`. It needs AVX-512 CPUs. The SLURM
  profile requests cascadelake nodes for its two rules for this reason.
- gemBS-rs `dev` at `baf6056` (`bs_call` 3.1.2, `mextr` 2.1.0, `read_filter`).
  `master` last changed in 2022. None of the three is packaged on conda.

Any change to a manifest, recipe, or lock reruns every rule.

## Outputs

Per sample, under `results/<run_id>/<sample>/`:

| Files | Step |
|---|---|
| `<s>.NNN.json` | gem-mapper report per part (`split/` FASTQs and part BAMs are temporary) |
| `<s>.bam`, `.bam.csi`, `.bam.md5` | Final alignments |
| `<s>.json` | Mapping QC: gem-mapper reports, `samtools stats` and average coverage |
| `<s>_coverage_{neg,pos}.bw` | Strand-specific coverage (the bedGraphs are temporary) |
| `<s>_<pool>.json` | Per-pool `bs_call` reports |
| `<s>.bcf`, `.bcf.csi`, `.bcf.md5` | Whole-genome methylation calls |
| `<s>_{cpg,non_cpg}.txt.gz`, `<s>_{chg,chh,cpg}.bed.gz`, `<s>_{neg,pos}.bw` | `mextr` outputs |
| `<s>_{cpg,chg,chh}_{pos,neg}.bw` | Methylation level per context and strand |
| `<s>_methylc.bed.gz` (+ `.tbi`) | Combined CG/CHG/CHH methylC file |

`results/refs/` holds `reference.fasta` (+ `.fai`), `chromsizes.tsv`, the pool
BEDs and their order in `pools/`,
`reference.contig_md5` (LN, M5, AS and SP per contig, as in gemBS) and the GEM
index `reference.gem`.

## Pipeline steps

1. **indexing**: decompress the reference, index it, write chromsizes and the contig metadata.
2. **gem_indexer**: build the GEM bisulfite index.
3. **split**: split each FASTQ pair with `seqkit split2`.
4. **mapping**: map each part with `gem-mapper`, clean read names and add the
   contig metadata with gemBS's `read_filter`, and sort, as gemBS does.
5. **merge**: merge the sorted parts into one read group and checksum the BAM.
6. **bam_qc**: write the mapping QC JSON with `samtools stats` and average coverage.
7. **get_coverage**: strand coverage with `bedtools genomecov` and `bedGraphToBigWig`.
8. **contig_pools**: group the contigs into `bs_call` jobs as gemBS does, once
   per reference. A contig of at least `calling.contig_pool_limit` bases is its
   own job. The rest are packed into `Pool@N` jobs of about that size.
9. **calling**: `bs_call` per pool, reading the BAM in place.
10. **merge_chroms**: concatenate the pool BCFs in reference order with `bcftools concat`.
11. **gembs_extract**: `mextr` bedMethyl, CpG/non-CpG and bigWig outputs.
12. **extract_signal**: methylation-level bigWigs per context and strand.
13. **get_methylc**: combine the three contexts into one sorted, indexed methylC file.

## Known issues

These come from gemBS-rs and were found by reading its code. They are not
tested here.

- `mextr` indexes its text outputs with tabix's generic preset, which reads the
  0-based start column as 1-based. Region queries on those `.tbi` files are off
  by one base.
- `mextr` rebuilds its regions as `contig:start-end` strings from 0-based
  coordinates and parses them as 1-based, which may drop the last base of each
  contig.

## References

Heath, Simon. n.d. GemBS-Rs. Github. Accessed October 10, 2026. https://github.com/heathsc/gemBS-rs.

Marco-Sola, Santiago, Michael Sammeth, Roderic Guigó, and Paolo Ribeca. 2012. “The GEM Mapper: Fast, Accurate and Versatile Alignment by Filtration.” Nature Methods 9 (12): 1185–88. https://doi.org/10.1038/nmeth.2221.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (April): 33. https://doi.org/10.12688/f1000research.29032.2.

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
