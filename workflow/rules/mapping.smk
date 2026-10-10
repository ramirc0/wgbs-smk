rule split:
    input:
        r1=lambda w: FASTQ_1_OF[w.sample],
        r2=lambda w: FASTQ_2_OF[w.sample],
    output:
        # seqkit -f wipes this directory, so it holds nothing else.
        dir=directory(f"{OUTDIR}/{{sample}}/split"),
        fastqs=expand(
            f"{OUTDIR}/{{{{sample}}}}/split/MOHD_{{{{sample}}}}_{{read}}.part_{{part}}.fastq.gz",
            part=SPLIT,
            read=["R1", "R2"],
        ),
    log:
        f"{LOGDIR}/split/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/split/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        echo "$(date): Rule started"

        seqkit split2 \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            -O {output.dir:q} \
            -p 2 \
            -j {threads} \
            -f

        echo "$(date): Rule finished"
        """


rule mapping:
    input:
        gem_index=f"{REFDIR}/reference.gem",
        r1=f"{OUTDIR}/{{sample}}/split/MOHD_{{sample}}_R1.part_{{part}}.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/split/MOHD_{{sample}}_R2.part_{{part}}.fastq.gz",
    output:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.{{part}}.bam",
        json=f"{OUTDIR}/{{sample}}/{{sample}}.{{part}}.json",
    params:
        read_group=lambda w: read_group(w.sample),
    log:
        f"{LOGDIR}/mapping/{{sample}}.{{part}}.txt",
    benchmark:
        f"{BENCHDIR}/mapping/{{sample}}.{{part}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        echo "$(date): Rule started"

        gem-mapper \
            -t {threads} \
            -I {input.gem_index:q} \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            -p \
            --report-file={output.json:q} \
            -r {params.read_group:q} > {output.bam:q}

        echo "$(date): Rule finished"
        """


rule merge:
    input:
        bams=expand(f"{OUTDIR}/{{{{sample}}}}/{{{{sample}}}}.{{part}}.bam", part=SPLIT),
    output:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.merged.bam",
    log:
        f"{LOGDIR}/merge/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/merge/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        echo "$(date): Rule started"

        samtools merge \
            -@ {threads} \
            -o - \
            {input.bams:q} | \
        samtools sort \
            -@ {threads} \
            --write-index \
            -T /tmp/{wildcards.sample:q} \
            -o {output.bam:q} -

        echo "$(date): Rule finished"
        """


rule postprocess:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.merged.bam",
        jsons=expand(f"{OUTDIR}/{{{{sample}}}}/{{{{sample}}}}.{{part}}.json", part=SPLIT),
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        bam_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bam.csi",
        md5=f"{OUTDIR}/{{sample}}/{{sample}}.bam.md5",
    log:
        f"{LOGDIR}/postprocess/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/postprocess/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        echo "$(date): Rule started"

        python workflow/scripts/read_filter.py \
            -i {input.bam:q} \
            {input.chromsizes:q} | \
        samtools sort \
            -o {output.bam:q} \
            -T /tmp/{wildcards.sample:q} \
            --threads {threads} \
            --write-index -
        echo "$(date): gem-mapper completed"

        md5sum {output.bam:q} > {output.md5:q}
        echo "$(date): md5sum completed"

        python workflow/scripts/make_average_coverage.py \
            --bamfile {output.bam:q} \
            --chromsizes {input.chromsizes:q} \
            --threads {threads} \
            --gem_mapper_jsons {input.jsons:q}
        echo "$(date): average coverage completed"

        echo "$(date): Rule finished"
        """


rule get_coverage:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        coverage_bedgraphs=[
            f"{OUTDIR}/{{sample}}/{{sample}}_coverage_neg.bg",
            f"{OUTDIR}/{{sample}}/{{sample}}_coverage_pos.bg",
        ],
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
