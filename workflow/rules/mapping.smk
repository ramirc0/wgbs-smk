rule split:
    input:
        r1=lambda w: SHEET[w.sample]["fastq_1"],
        r2=lambda w: SHEET[w.sample]["fastq_2"],
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

        # The links fix the part names: seqkit keeps each input's name and extension.
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
        row = SHEET[wildcards.sample]
        return {"r1": row["fastq_1"], "r2": row["fastq_2"]}
    split = f"{OUTDIR}/{wildcards.sample}/split/{wildcards.sample}"
    return {
        "r1": f"{split}_R1.part_{wildcards.part}.fastq.gz",
        "r2": f"{split}_R2.part_{wildcards.part}.fastq.gz",
    }


# One part maps straight to the final BAM, as gemBS maps a whole dataset.
STEM = "{sample}" if len(SPLIT) == 1 else "{sample}.{part}"

# gem-mapper writes a report per job. minibwa writes none.
MAP_REPORTS = [f"{OUTDIR}/{{sample}}/{STEM}.map.json"] if ALIGNER == "gem3" else []

INDEX = {
    "gem3": [f"{REFDIR}/reference.gem"],
    "minibwa": [f"{REFDIR}/reference.l2b", f"{REFDIR}/reference.meth.mbw"],
}[ALIGNER]


rule:
    # Each aligner's jobs get their own rule name, so the profiles can size them.
    name:
        "mapping_minibwa" if ALIGNER == "minibwa" else "mapping"
    input:
        unpack(mapping_fastqs),
        index=INDEX,
        contig_md5=f"{REFDIR}/reference.contig_md5",
    output:
        bam=(temp if len(SPLIT) > 1 else str)(f"{OUTDIR}/{{sample}}/{STEM}.bam"),
        bam_csi=(temp if len(SPLIT) > 1 else str)(f"{OUTDIR}/{{sample}}/{STEM}.bam.csi"),
        json=MAP_REPORTS,
    params:
        aligner=ALIGNER,
        read_group=lambda w: read_group(w.sample),
        tmpdir=config["tmpdir"],
    log:
        f"{LOGDIR}/mapping/{STEM}.txt",
    benchmark:
        f"{BENCHDIR}/mapping/{STEM}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        mkdir -p {params.tmpdir:q}
        if [[ {params.aligner:q} == gem3 ]]; then
            gem-mapper \
                --threads {threads} \
                -I {input.index:q} \
                --i1 {input.r1:q} \
                --i2 {input.r2:q} \
                --paired-end-alignment \
                --bisulfite-conversion inferred-C2T-G2A \
                --report-file {output.json:q} \
                --sam-read-group-header {params.read_group:q}
        else
            # bs_call counts a read without a strand tag as unconverted, and
            # minibwa writes none. In a directional library, R1 forward and R2
            # reverse are C-to-T: gem-mapper's XB:A:C.
            minibwa map \
                --meth \
                -t {threads} \
                -R {params.read_group:q} \
                $(dirname {input.index[0]:q})/reference \
                {input.r1:q} \
                {input.r2:q} \
            | awk '/^@/ {{print; next}} {{print $0 "\tXB:A:" (int($2 / 16) % 2 == int($2 / 128) % 2 ? "C" : "G")}}'
        fi \
        | read_filter {input.contig_md5:q} \
        | samtools sort \
            -@ {threads} \
            -T {params.tmpdir:q}/$(basename {output.bam:q}) \
            --write-index \
            -o {output.bam:q}##idx##{output.bam_csi:q} \
            -
        """


if len(SPLIT) > 1:

    rule merge:
        input:
            bams=expand(f"{OUTDIR}/{{sample}}/{STEM}.bam", part=SPLIT, allow_missing=True),
        output:
            bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
            bam_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bam.csi",
        log:
            f"{LOGDIR}/merge/{{sample}}.txt",
        benchmark:
            f"{BENCHDIR}/merge/{{sample}}.tsv",
        software:
            SOFTWARE_ENV
        shell:
            r"""
            exec &> >(tee {log:q})

            # -c and -p keep one @RG and @PG per dataset: every part has the same IDs.
            samtools merge \
                -c \
                -p \
                -@ {threads} \
                --write-index \
                -o {output.bam:q}##idx##{output.bam_csi:q} \
                {input.bams:q}
            """


rule bam_md5:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
    output:
        md5=f"{OUTDIR}/{{sample}}/{{sample}}.bam.md5",
    log:
        f"{LOGDIR}/bam_md5/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bam_md5/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        md5sum {input.bam:q} > {output.md5:q}
        """


rule bam_qc:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        jsons=expand(MAP_REPORTS, part=SPLIT, allow_missing=True),
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
        coverage_bigwigs=[
            f"{OUTDIR}/{{sample}}/{{sample}}_coverage_neg.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_coverage_pos.bw",
        ],
    params:
        tmpdir=config["tmpdir"],
    log:
        f"{LOGDIR}/get_coverage/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/get_coverage/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # One genomecov per contig and strand, largest contigs first so the
        # last jobs are short. A contig's two strands run side by side, so the
        # second reads its BAM region from the page cache.
        mkdir -p {params.tmpdir:q}
        STAGE=$(mktemp -d -p {params.tmpdir:q})
        mkdir "$STAGE"/neg "$STAGE"/pos
        parallel \
            -j {threads} \
            "set -o pipefail; samtools view -u {input.bam:q} {{1}} | bedtools genomecov -ibam stdin -bg -strand {{2}} > $STAGE/{{3}}/{{1}}.bg" \
            :::: <(sort -k2,2nr {input.chromsizes:q} | cut -f1) \
            ::: - + \
            :::+ neg pos

        # bigtools takes contigs in sort -k1,1 order. --inmemory keeps it off
        # disk, where a failed write shows only as SendError::Disconnected.
        cut -f1 {input.chromsizes:q} | LC_ALL=C sort > "$STAGE"/order.txt
        parallel \
            -j 2 \
            "set -o pipefail; sed 's|.*|$STAGE/{{1}}/&.bg|' $STAGE/order.txt | xargs -d '\n' cat | bigtools bedgraphtobigwig --inmemory -t $(( ({threads} + 1) / 2 )) - {input.chromsizes:q} {{2}}" \
            ::: neg pos \
            :::+ {output.coverage_bigwigs:q}
        rm -r "$STAGE"
        """
