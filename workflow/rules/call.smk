rule generate_chrome_ctgs:
    input:
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        ctg_beds=temp(
            expand(
                f"{OUTDIR}/{{{{sample}}}}/{{{{sample}}}}_{{contig}}_ctgs.bed",
                contig=CONTIGS,
            )
        ),
    params:
        stage=lambda w: f"{config['tmpdir']}/{OUTDIR}/calls/{w.sample}",
    log:
        f"{LOGDIR}/generate_chrome_ctgs/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/generate_chrome_ctgs/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        OUTPUT_DIR={params.stage:q}
        mkdir -p "$OUTPUT_DIR"
        cp {input.chromsizes:q} $OUTPUT_DIR

        python workflow/scripts/make_ctg_beds.py \
            --chromsizes $OUTPUT_DIR/$(basename {input.chromsizes:q}) \
            --output_dir $OUTPUT_DIR \
            --barcode {wildcards.sample:q}

        rm -rf $OUTPUT_DIR/$(basename {input.chromsizes:q})
        cp -t $(dirname {output.ctg_beds[0]:q}) "$OUTPUT_DIR"/{wildcards.sample:q}_*_ctgs.bed
        rm -rf "$OUTPUT_DIR"

        echo "Generated contig bed files!"
        """


rule calling:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.bam",
        bam_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bam.csi",
        ctg_bed=f"{OUTDIR}/{{sample}}/{{sample}}_{{contig}}_ctgs.bed",
        fasta=f"{REFDIR}/reference.fasta",
        fasta_fai=f"{REFDIR}/reference.fasta.fai",
        chromsizes=f"{REFDIR}/chromsizes.tsv",
    output:
        bcf=temp(f"{OUTDIR}/{{sample}}/{{sample}}_{{contig}}.bcf"),
        report=f"{OUTDIR}/{{sample}}/{{sample}}_{{contig}}.json",
    params:
        stage=lambda w: f"{config['tmpdir']}/{OUTDIR}/calls/{w.sample}",
        out_type=lambda w: config["calling"]["output_type"],
        conversion=lambda w: config["calling"]["conversion"],
        ltrim=lambda w: config["calling"]["left_trim"],
        rtrim=lambda w: config["calling"]["right_trim"],
        ref_bias=lambda w: config["calling"]["reference_bias"],
        mapq=lambda w: config["calling"]["mapq_threshold"],
        qual=lambda w: config["calling"]["qual_threshold"],
    log:
        f"{LOGDIR}/calling/{{sample}}_{{contig}}.txt",
    benchmark:
        f"{BENCHDIR}/calling/{{sample}}_{{contig}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        OUTPUT_DIR={params.stage:q}/{wildcards.contig:q}
        CTGS_DIR={params.stage:q}
        mkdir -p "$OUTPUT_DIR"
        mkdir -p "$CTGS_DIR"
        cp {input.bam:q} {input.bam_csi:q} {input.chromsizes:q} {input.fasta:q} {input.fasta_fai:q} $OUTPUT_DIR
        cp {input.ctg_bed:q} $CTGS_DIR

        bs_call \
            --loglevel info \
            --output $OUTPUT_DIR/$(basename {output.bcf:q}) \
            --output-type {params.out_type:q} \
            --reference $OUTPUT_DIR/$(basename {input.fasta:q}) \
            --sample {wildcards.sample:q} \
            --contig-include $OUTPUT_DIR/$(basename {input.chromsizes:q}) \
            --report-file $OUTPUT_DIR/$(basename {output.report:q}) \
            --contig-bed $CTGS_DIR/$(basename {input.ctg_bed:q}) \
            --threads {threads} \
            --conversion {params.conversion:q} \
            --left-trim {params.ltrim:q} \
            --right-trim {params.rtrim:q} \
            --reference-bias {params.ref_bias:q} \
            --mapq-threshold {params.mapq:q} \
            --bq-threshold {params.qual:q} \
            $OUTPUT_DIR/$(basename {input.bam:q})

        cp $OUTPUT_DIR/$(basename {output.bcf:q}) $OUTPUT_DIR/$(basename {output.report:q}) $(dirname {output.bcf:q})
        rm -rf $OUTPUT_DIR/$(basename {input.fasta:q}) $OUTPUT_DIR/$(basename {input.fasta:q}) $OUTPUT_DIR/$(basename {input.bam:q}) $OUTPUT_DIR/$(basename {input.bam_csi:q})
        rm -rf $CTGS_DIR/$(basename {input.ctg_bed:q})

        echo "Generated BCF file for {wildcards.contig:q}"
        """


rule merge_chroms:
    input:
        chrom_bcfs=expand(
            f"{OUTDIR}/{{{{sample}}}}/{{{{sample}}}}_{{contig}}.bcf", contig=CONTIGS
        ),
    output:
        bcf=f"{OUTDIR}/{{sample}}/{{sample}}.bcf",
        bcf_csi=f"{OUTDIR}/{{sample}}/{{sample}}.bcf.csi",
        bcf_md5=f"{OUTDIR}/{{sample}}/{{sample}}.bcf.md5",
    params:
        stage=lambda w: f"{config['tmpdir']}/{OUTDIR}/calls/{w.sample}",
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

        OUTPUT_DIR={params.stage:q}
        mkdir -p "$OUTPUT_DIR"
        cp {input.chrom_bcfs:q} $OUTPUT_DIR

        # Contig order in the merged BCF is the shell's glob order, not CONTIGS order.
        bcftools concat \
            --output $OUTPUT_DIR/$(basename {output.bcf:q}) \
            --output-type {params.out_type:q} \
            --naive \
            --threads {threads} \
            $OUTPUT_DIR/{wildcards.sample:q}_*.bcf

        bcftools index --threads {threads} $OUTPUT_DIR/$(basename {output.bcf:q})

        md5sum $OUTPUT_DIR/$(basename {output.bcf:q}) > $OUTPUT_DIR/$(basename {output.bcf_md5:q})

        cp -t $(dirname {output.bcf:q}) $OUTPUT_DIR/$(basename {output.bcf:q}) $OUTPUT_DIR/$(basename {output.bcf_csi:q}) $OUTPUT_DIR/$(basename {output.bcf_md5:q})
        rm -rf "$OUTPUT_DIR"

        echo "Done!"
        """
