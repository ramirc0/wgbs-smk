rule split:
    input:
        r1=lambda w: FASTQ_1_OF[w.sample],
        r2=lambda w: FASTQ_2_OF[w.sample],
    output:
        # seqkit -f wipes this directory, so it holds nothing else.
        dir=directory(f"{OUTDIR}/{{sample}}/split"),
        fastqs=temp(
            expand(
                f"{OUTDIR}/{{{{sample}}}}/split/{{{{sample}}}}_{{read}}.part_{{part}}.fastq.gz",
                part=SPLIT,
                read=["R1", "R2"],
            )
        ),
    params:
        parts=len(SPLIT),
        tmpdir=config["tmpdir"],
    log:
        f"{LOGDIR}/split/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/split/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # seqkit names the parts after its inputs. --by-part-prefix is a no-op in 2.11.
        mkdir -p {params.tmpdir:q}
        STAGE=$(mktemp -d -p {params.tmpdir:q})
        ln -s "$(realpath {input.r1:q})" "$STAGE"/{wildcards.sample:q}_R1.fastq.gz
        ln -s "$(realpath {input.r2:q})" "$STAGE"/{wildcards.sample:q}_R2.fastq.gz

        seqkit split2 \
            -1 "$STAGE"/{wildcards.sample:q}_R1.fastq.gz \
            -2 "$STAGE"/{wildcards.sample:q}_R2.fastq.gz \
            -O {output.dir:q} \
            -p {params.parts} \
            -j {threads} \
            -f
        rm -r "$STAGE"
        """


def mapping_fastqs(wildcards):
    """Return the FASTQ pair one mapping job reads.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Job wildcards with `sample` and `part`.

    Returns
    -------
    dict of str to str
        `r1` and `r2`: the sheet's FASTQs when `fastq_split` is 1, as gemBS maps
        a dataset whole, else the split parts.
    """
    if len(SPLIT) == 1:
        return {"r1": FASTQ_1_OF[wildcards.sample], "r2": FASTQ_2_OF[wildcards.sample]}
    split = f"{OUTDIR}/{wildcards.sample}/split/{wildcards.sample}"
    return {
        "r1": f"{split}_R1.part_{wildcards.part}.fastq.gz",
        "r2": f"{split}_R2.part_{wildcards.part}.fastq.gz",
    }


rule mapping:
    input:
        unpack(mapping_fastqs),
        gem_index=f"{REFDIR}/reference.gem",
        contig_md5=f"{REFDIR}/reference.contig_md5",
    output:
        bam=temp(f"{OUTDIR}/{{sample}}/{{sample}}.{{part}}.bam"),
        bam_csi=temp(f"{OUTDIR}/{{sample}}/{{sample}}.{{part}}.bam.csi"),
        json=f"{OUTDIR}/{{sample}}/{{sample}}.{{part}}.json",
    params:
        read_group=lambda w: read_group(w.sample),
        sort_prefix=lambda w: f"{config['tmpdir']}/{w.sample}.{w.part}",
    log:
        f"{LOGDIR}/mapping/{{sample}}.{{part}}.txt",
    benchmark:
        f"{BENCHDIR}/mapping/{{sample}}.{{part}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        mkdir -p $(dirname {params.sort_prefix:q})
        gem-mapper \
            -t {threads} \
            -I {input.gem_index:q} \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            -p \
            --report-file={output.json:q} \
            -r {params.read_group:q} \
        | read_filter {input.contig_md5:q} \
        | samtools sort \
            -@ {threads} \
            -T {params.sort_prefix:q} \
            --write-index \
            -o {output.bam:q}##idx##{output.bam_csi:q} \
            -
        """


rule merge:
    input:
        bams=expand(f"{OUTDIR}/{{{{sample}}}}/{{{{sample}}}}.{{part}}.bam", part=SPLIT),
        csis=expand(f"{OUTDIR}/{{{{sample}}}}/{{{{sample}}}}.{{part}}.bam.csi", part=SPLIT),
    output:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        bam_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bam.csi",
        md5=f"{OUTDIR}/{{sample}}/{{sample}}.bam.md5",
    params:
        parts=len(SPLIT),
    log:
        f"{LOGDIR}/merge/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/merge/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        if [[ {params.parts} -eq 1 ]]; then
            ln -f {input.bams:q} {output.bam:q}
            ln -f {input.csis:q} {output.bam_csi:q}
        else
            # -c and -p keep one @RG and @PG per dataset: every part has the same IDs.
            samtools merge \
                -c \
                -p \
                -@ {threads} \
                --write-index \
                -o {output.bam:q}##idx##{output.bam_csi:q} \
                {input.bams:q}
        fi
        md5sum {output.bam:q} > {output.md5:q}
        """


rule bam_qc:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        jsons=expand(f"{OUTDIR}/{{{{sample}}}}/{{{{sample}}}}.{{part}}.json", part=SPLIT),
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        json=f"{OUTDIR}/{{sample}}/{{sample}}.json",
    log:
        f"{LOGDIR}/bam_qc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bam_qc/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/make_average_coverage.py \
            --bamfile {input.bam:q} \
            --chromsizes {input.chromsizes:q} \
            --threads {threads} \
            --gem_mapper_jsons {input.jsons:q} \
            --output {output.json:q}
        """


rule get_coverage:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        coverage_bedgraphs=temp(
            [
                f"{OUTDIR}/{{sample}}/{{sample}}_coverage_neg.bg",
                f"{OUTDIR}/{{sample}}/{{sample}}_coverage_pos.bg",
            ]
        ),
        coverage_bigwigs=[
            f"{OUTDIR}/{{sample}}/{{sample}}_coverage_neg.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_coverage_pos.bw",
        ],
    log:
        f"{LOGDIR}/get_coverage/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/get_coverage/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        echo "$(date): Rule started"

        echo "$(date --iso=seconds): Running bedtools genomecov"
        parallel \
            -j 2 \
            bedtools genomecov -ibam {input.bam:q} -bg -strand {{1}} '>' {{2}} \
            ::: - + \
            :::+ {output.coverage_bedgraphs:q}

        echo "$(date --iso=seconds): Running bedGraphToBigWig"
        parallel \
            -j 2 \
            bedGraphToBigWig {{1}} {input.chromsizes:q} {{2}} \
            ::: {output.coverage_bedgraphs:q} \
            :::+ {output.coverage_bigwigs:q}

        echo "$(date): Rule finished"
        """
