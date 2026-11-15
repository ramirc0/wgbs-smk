# The run's QC report: every sample's QC in one table, graded on ENCODE's WGBS
# portal audits (workflow/report/qc.yaml), with interactive figures.

import subprocess

REPORT_DIR = f"{OUTDIR}/report"
# The header names the pipeline commit, so a new commit renders the page again.
PIPELINE_VERSION = subprocess.run(
    ["git", "-C", workflow.basedir, "describe", "--always", "--dirty"],
    capture_output=True,
    text=True,
).stdout.strip()


def call_reports(wildcards):
    """Return every sample's per-pool bs_call reports.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Job wildcards, unused: the report covers every sample.

    Returns
    -------
    list of str
        One JSON per sample and contig pool.
    """
    return [f"{prefix(s, 'call')}_{pool}.json" for s in SAMPLES for pool in pool_names()]


rule gather_qc:
    input:
        mapping_qc=expand(f"{OUTDIR}/{{sample}}/align/{{sample}}_qc.json", sample=SAMPLES),
        call_reports=call_reports,
        pools=f"{REFDIR}/pools",
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        f"{REPORT_DIR}/qc_metrics.tsv",
    params:
        run=OUTDIR,
        samples=SAMPLES,
    log:
        f"{LOGDIR}/gather_qc/all.txt",
    benchmark:
        f"{BENCHDIR}/gather_qc/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/gather_qc.py \
            --run-dir {params.run:q} \
            --samples {params.samples:q} \
            --pools {input.pools:q} \
            --chromsizes {input.chromsizes:q} \
            --threads {threads:q} \
            --output {output:q}

        sync {output:q}
        """


rule qc_figures:
    input:
        mapping_qc=expand(f"{OUTDIR}/{{sample}}/align/{{sample}}_qc.json", sample=SAMPLES),
        call_reports=call_reports,
        pools=f"{REFDIR}/pools",
    output:
        expand(f"{REPORT_DIR}/samples/{{sample}}.js", sample=SAMPLES),
    params:
        run=OUTDIR,
        report_dir=REPORT_DIR,
        samples=SAMPLES,
    log:
        f"{LOGDIR}/qc_figures/all.txt",
    benchmark:
        f"{BENCHDIR}/qc_figures/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/qc_figures.py \
            --run-dir {params.run:q} \
            --report-dir {params.report_dir:q} \
            --samples {params.samples:q} \
            --pools {input.pools:q} \
            --threads {threads:q}

        sync {output:q}
        """


rule qc_report:
    input:
        spec=f"{workflow.basedir}/report/qc.yaml",
        metrics=f"{REPORT_DIR}/qc_metrics.tsv",
        figures=expand(f"{REPORT_DIR}/samples/{{sample}}.js", sample=SAMPLES),
    output:
        f"{REPORT_DIR}/qc_report.html",
    params:
        run=OUTDIR,
        run_id=RUN_ID,
        version=f"wgbs-smk {PIPELINE_VERSION}",
    log:
        f"{LOGDIR}/qc_report/all.txt",
    benchmark:
        f"{BENCHDIR}/qc_report/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        qc-report \
            --spec {input.spec:q} \
            --metrics {input.metrics:q} \
            --run-dir {params.run:q} \
            --run-id {params.run_id:q} \
            --version {params.version:q} \
            --output {output:q}

        sync {output:q}
        """
