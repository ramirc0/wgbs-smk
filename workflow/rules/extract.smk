rule gembs_extract:
    input:
        bcf=f"{OUTDIR}/{{sample}}/{{sample}}.bcf",
        bcf_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bcf.csi",
    output:
        cpgfile=f"{OUTDIR}/{{sample}}/{{sample}}_cpg.txt.gz",
        noncpgfile=f"{OUTDIR}/{{sample}}/{{sample}}_non_cpg.txt.gz",
        bed_gz=[
            f"{OUTDIR}/{{sample}}/{{sample}}_chg.bed.gz",
            f"{OUTDIR}/{{sample}}/{{sample}}_chh.bed.gz",
            f"{OUTDIR}/{{sample}}/{{sample}}_cpg.bed.gz",
        ],
        bigwigs=[
            f"{OUTDIR}/{{sample}}/{{sample}}_neg.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_pos.bw",
        ],
    params:
        # Separate from the calls staging dir: this rule copies back everything it finds.
        stage=lambda w: f"{config['tmpdir']}/{OUTDIR}/extract/{w.sample}",
        flags=lambda w: extract_flags(),
    log:
        f"{LOGDIR}/gembs_extract/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/gembs_extract/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        OUTPUT_DIR={params.stage:q}
        mkdir -p "$OUTPUT_DIR"

        cp {input.bcf:q} {input.bcf_csi:q} $OUTPUT_DIR

        mextr \
            --loglevel info \
            --compress \
            --md5 \
            --cpgfile $OUTPUT_DIR/$(basename {output.cpgfile:q}) \
            --noncpgfile $OUTPUT_DIR/$(basename {output.noncpgfile:q}) \
            --bed-methyl $OUTPUT_DIR/{wildcards.sample:q} \
            --tabix \
            --threads {threads} \
            {params.flags:q} \
            $OUTPUT_DIR/$(basename {input.bcf:q})

        rm -rf $OUTPUT_DIR/$(basename {input.bcf:q}) $OUTPUT_DIR/$(basename {input.bcf_csi:q})
        cp -t $(dirname {output.cpgfile:q}) "$OUTPUT_DIR"/*
        rm -rf "$OUTPUT_DIR"
        """


rule extract_signal:
    input:
        chromsizes=f"{REFDIR}/chromsizes.tsv",
        bed_gz=[
            f"{OUTDIR}/{{sample}}/{{sample}}_chg.bed.gz",
            f"{OUTDIR}/{{sample}}/{{sample}}_chh.bed.gz",
            f"{OUTDIR}/{{sample}}/{{sample}}_cpg.bed.gz",
        ],
    output:
        bigwigs=[
            f"{OUTDIR}/{{sample}}/{{sample}}_cpg_neg.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_cpg_pos.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_chg_pos.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_chg_neg.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_chh_pos.bw",
            f"{OUTDIR}/{{sample}}/{{sample}}_chh_neg.bw",
        ],
    log:
        f"{LOGDIR}/extract_signal/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/extract_signal/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/extract_signal.py \
            --bed-gz {input.bed_gz:q} \
            --chromsizes {input.chromsizes:q} \
            --output-dir $(dirname {output.bigwigs[0]:q}) \
            --barcode {wildcards.sample:q}
        """


rule get_methylc:
    input:
        chh_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_chh.bed.gz",
        chg_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_chg.bed.gz",
        cpg_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_cpg.bed.gz",
    output:
        methylc_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_methylc.bed.gz",
    log:
        f"{LOGDIR}/get_methylc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/get_methylc/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/get_methylc.py \
            --barcode {wildcards.sample:q} \
            --chh {input.chh_bed_gz:q} \
            --chg {input.chg_bed_gz:q} \
            --cpg {input.cpg_bed_gz:q} \
            --threads {threads} \
            --outdir $(dirname {output.methylc_bed_gz:q})
        """
