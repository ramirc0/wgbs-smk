# WGBS Snakemake Pipeline

A Snakemake workflow for whole genome bisulfite sequencing (WGBS). It maps
paired-end reads with [GEM3](https://github.com/smarco/gem3-mapper), calls
methylation with [gemBS-rs](https://github.com/heathsc/gemBS-rs) `bs_call`, and
extracts bedMethyl, methylC and bigWig tracks with `mextr`. It follows the
`snakemake-template` layout and conventions.

## Layout

```
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per sample: sample_id, dataset, fastq_1, fastq_2
workflow/
  Snakefile              includes the rules + rule all
  rules/
    common.smk           config, sample sheet, path constants, helpers
    indexing.smk         reference FASTA, chromsizes, SAM dict, GEM index
    mapping.smk          split, map, merge, filter, strand coverage
    call.smk             per-contig bs_call, merged BCF
    extract.smk          mextr, strand bigWigs, methylC
  scripts/               argparse scripts, each with build_parser() + main()
  envs/
    pixi.toml            rule env (pixi workspace)
    recipes/             gem3-mapper and gemBS-rs, built from source
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + lint task
pyproject.toml           ruff + numpydoc config
```

## Idioms

- **Outputs grouped by `run_id`** under `results/<run_id>/<sample>/`. Logs and
  benchmarks mirror it (`logs/<run_id>/`, `benchmarks/<run_id>/`). The
  reference artifacts are run-independent and go in `results/refs/`.
- **Sample sheet** is a TSV read once in `common.smk`. Columns are looked up
  with `_column()` (empty cell -> `None`).
- **Config -> flags**: helpers such as `extract_flags()` and `dict_flags()` turn
  `config` into flag lists. `read_group()` builds the gem-mapper `@RG` line.
- **Every rule** carries `log:`, `benchmark:`, and `software:`.
- **Staging**: the heavy rules copy their inputs to `tmpdir`, run there, and
  copy the results back.

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
```

Run from the repository root. Local runs suit small data only.

The first job builds the rule env, which compiles gem3-mapper and gemBS-rs.
Build it once before launching many jobs:

```bash
pixi install --manifest-path workflow/envs/pixi.toml
```

Tune `set-resources` in the profiles from the `benchmarks/<run_id>/` TSVs after
a first run.

### Sample sheet

One row per sample. Lines starting with `#` are ignored.

- `sample_id`: unique ID, letters and digits only.
- `dataset`: dataset or run name. It fills the read group `ID` and `PU`.
- `fastq_1`, `fastq_2`: paired FASTQ paths. Their basenames MUST be
  `MOHD_<sample_id>_R1.fastq.gz` and `MOHD_<sample_id>_R2.fastq.gz` (see
  known issues).

### Configuration

- `samples`, `outdir`, `run_id`, `pixi_env`: as in the template.
- `tmpdir`: node-local scratch for staging. On systems where `/tmp` is mounted
  `noexec`, pick another path. The SLURM profile also sets `shadow-prefix: /tmp`.
- `references.reference`: reference FASTA, plain or gzipped.
- `assembly`, `species` (optional): `AS` and `SP` tags of the SAM sequence
  dictionary. `species` MUST NOT contain whitespace.
- `index.sampling_rate`: gem-indexer text sampling rate. Lower values give a
  larger, faster index.
- `fastq_split`: number of parts per FASTQ pair. Parts are mapped in parallel
  and merged.
- `calling.*`: `bs_call` output type, MAPQ and base-quality thresholds,
  reference bias, read trimming, and conversion rates
  (`under,over`, comma-delimited).
- `extract.*`: `mextr` strand-specific bedMethyl and bigWig modes, and the
  Phred threshold.

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main, the pixi
software-deployment plugin and the SLURM executor plugins. The rule env is a
second pixi workspace, `workflow/envs/`. Its `pixi.lock` pins every package.
pixi MUST be on `PATH`.

Most tools come from conda-forge and bioconda. Two are built from source by
`pixi-build`, at the commits the earlier Docker image used:

- gem3-mapper `master` (27 commits past the 3.6.1 release), compiled with
  `-O3 -march=x86-64-v4`. It needs AVX-512 CPUs. The SLURM profile requests
  cascadelake nodes for this reason.
- gemBS-rs `master` (`bs_call`, `mextr`). Neither is packaged on conda.

Any change to a manifest, recipe, or lock reruns every rule.

## Outputs

Per sample, under `results/<run_id>/<sample>/`:

| Files | Step |
|---|---|
| `split/MOHD_<s>_R{1,2}.part_NNN.fastq.gz` | FASTQ parts for parallel mapping |
| `<s>.NNN.bam`, `<s>.NNN.json` | Per-part alignments (SAM text) and gem-mapper reports |
| `<s>.merged.bam` | Merged, sorted alignments |
| `<s>.bam`, `.bam.csi`, `.bam.md5`, `<s>.json` | Final alignments and the mapping QC JSON with average coverage |
| `<s>_coverage_{neg,pos}.{bg,bw}` | Strand-specific coverage |
| `<s>_<contig>.json` | Per-contig `bs_call` reports |
| `<s>.bcf`, `.bcf.csi`, `.bcf.md5` | Whole-genome methylation calls |
| `<s>_{cpg,non_cpg}.txt.gz`, `<s>_{chg,chh,cpg}.bed.gz`, `<s>_{neg,pos}.bw` | `mextr` outputs |
| `<s>_{cpg,chg,chh}_{pos,neg}.bw` | Methylation level per context and strand |
| `<s>_methylc.bed.gz` (+ `.tbi`) | Combined CG/CHG/CHH methylC file |

`results/refs/` holds `reference.fasta` (+ `.fai`), `chromsizes.tsv`,
`sam_dict.tsv` and the GEM index `reference.gem`.

## Pipeline steps

1. **indexing**: decompress the reference, index it, write chromsizes and the SAM dictionary.
2. **gem_indexer**: build the GEM bisulfite index.
3. **split**: split each FASTQ pair with `seqkit split2`.
4. **mapping**: map each part with `gem-mapper`.
5. **merge**: merge and sort the parts.
6. **postprocess**: clean read names, sort, checksum, and add average coverage to the QC JSON.
7. **get_coverage**: strand coverage with `bedtools genomecov` and `bedGraphToBigWig`.
8. **generate_chrome_ctgs**: one contig BED per main chromosome, plus a `Pool` of the rest.
9. **calling**: `bs_call` per contig.
10. **merge_chroms**: concatenate the contig BCFs with `bcftools concat`.
11. **gembs_extract**: `mextr` bedMethyl, CpG/non-CpG and bigWig outputs.
12. **extract_signal**: methylation-level bigWigs per context and strand.
13. **get_methylc**: combine the three contexts into one sorted, indexed methylC file.

## Known issues

These behaviors carry over unchanged from the pre-template pipeline (git tag
`pre-template`).

- `split` always passes `-p 2` to `seqkit split2`. Any `fastq_split` other
  than 2 fails with missing outputs.
- `split` outputs assume FASTQ basenames `MOHD_<sample_id>_R{1,2}.fastq.gz`.
  Other names fail with missing outputs.
- `indexing` fails if the reference's basename is `reference.fasta`: `mv` onto
  itself errors.
- `CONTIGS` hardcodes human chromosomes (chr1-22, X, Y, Pool). A reference
  missing any of them fails at `generate_chrome_ctgs`.
- Contig BEDs depend only on chromsizes but are rebuilt per sample.
- The contig order in the merged BCF follows the shell's glob order, which
  depends on the locale.
- `calling` removes the staged FASTA twice and leaves the FASTA index,
  chromsizes and its outputs in `tmpdir`.
- Several files are written but not declared as outputs: the merged BAM's
  `.csi`, `<s>.json`, the methylC `.tbi`, and the `mextr` `.md5`, `.tbi` and
  `.bb` files.
- `extract_signal` writes an empty file in place of a bigWig for a strand with
  no sites.
- Split FASTQs and per-part BAMs are kept, which doubles the input size on disk.

## References

Heath, Simon. n.d. GemBS-Rs. Github. Accessed October 10, 2026. https://github.com/heathsc/gemBS-rs.

Marco-Sola, Santiago, Michael Sammeth, Roderic Guigó, and Paolo Ribeca. 2012. “The GEM Mapper: Fast, Accurate and Versatile Alignment by Filtration.” Nature Methods 9 (12): 1185–88. https://doi.org/10.1038/nmeth.2221.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (April): 33. https://doi.org/10.12688/f1000research.29032.2.

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
