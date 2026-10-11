checkpoint contig_pools:
    input:
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        pools=directory(f"{REFDIR}/pools"),
    params:
        limit=lambda w: config["calling"].get("contig_pool_limit", 25_000_000),
    log:
        f"{LOGDIR}/contig_pools.txt",
    benchmark:
        f"{BENCHDIR}/contig_pools.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/contig_pools.py \
            --chromsizes {input.chromsizes:q} \
            --limit {params.limit:q} \
            --output-dir {output.pools:q}
        """


rule calling:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        bam_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bam.csi",
        pools=f"{REFDIR}/pools",
        fasta=f"{REFDIR}/reference.fasta",
        fasta_fai=f"{REFDIR}/reference.fasta.fai",
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        bcf=temp(f"{OUTDIR}/{{sample}}/{{sample}}_{{pool}}.bcf"),
        report=f"{OUTDIR}/{{sample}}/{{sample}}_{{pool}}.json",
    params:
        out_type=lambda w: config["calling"]["output_type"],
        flags=lambda w: call_flags(),
    log:
        f"{LOGDIR}/calling/{{sample}}_{{pool}}.txt",
    benchmark:
        f"{BENCHDIR}/calling/{{sample}}_{{pool}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bs_call \
            --loglevel info \
            --output {output.bcf:q} \
            --output-type {params.out_type:q} \
            --reference {input.fasta:q} \
            --sample {wildcards.sample:q} \
            --contig-include {input.chromsizes:q} \
            --report-file {output.report:q} \
            --contig-bed {input.pools:q}/{wildcards.pool:q}.bed \
            --threads {threads} \
            {params.flags:q} \
            {input.bam:q}
        """


def pool_bcfs(wildcards):
    """Return a sample's per-pool BCFs in reference order.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Job wildcards with `sample`.

    Returns
    -------
    list of str
        One BCF per pool, in the order of `pools.txt`.
    """
    pools = Path(checkpoints.contig_pools.get().output.pools, "pools.txt")
    return [f"{prefix(wildcards.sample)}_{pool}.bcf" for pool in pools.read_text().split()]


rule merge_chroms:
    input:
        bcfs=pool_bcfs,
    output:
        bcf=f"{OUTDIR}/{{sample}}/{{sample}}.bcf",
        bcf_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bcf.csi",
        bcf_md5=f"{OUTDIR}/{{sample}}/{{sample}}.bcf.md5",
    params:
        out_type=lambda w: config["calling"]["output_type"],
    log:
        f"{LOGDIR}/merge_chroms/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/merge_chroms/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools concat \
            --output {output.bcf:q} \
            --output-type {params.out_type:q} \
            --naive \
            --threads {threads} \
            {input.bcfs:q}
        bcftools index --threads {threads} {output.bcf:q}
        md5sum {output.bcf:q} > {output.bcf_md5:q}
        """
