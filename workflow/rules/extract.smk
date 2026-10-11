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
        # mextr writes a .md5 for every output and a .tbi for every text one.
        sidecars=[
            f"{OUTDIR}/{{sample}}/{{sample}}_{name}{ext}"
            for name, exts in [
                ("cpg.txt.gz", (".tbi", ".md5")),
                ("non_cpg.txt.gz", (".tbi", ".md5")),
                ("chg.bed.gz", (".tbi", ".md5")),
                ("chh.bed.gz", (".tbi", ".md5")),
                ("cpg.bed.gz", (".tbi", ".md5")),
                ("chg.bb", ("", ".md5")),
                ("chh.bb", ("", ".md5")),
                ("cpg.bb", ("", ".md5")),
                ("neg.bw", (".md5",)),
                ("pos.bw", (".md5",)),
            ]
            for ext in exts
        ],
    params:
        flags=lambda w: extract_flags(),
        tmpdir=config["tmpdir"],
    log:
        f"{LOGDIR}/gembs_extract/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/gembs_extract/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # mextr's tabix thread reads each output back while it is written. On
        # CephFS those reads come back short, so mextr writes to local disk.
        mkdir -p {params.tmpdir:q}
        STAGE=$(mktemp -d -p {params.tmpdir:q})
        mextr \
            --loglevel info \
            --compress \
            --md5 \
            --cpgfile "$STAGE"/$(basename {output.cpgfile:q}) \
            --noncpgfile "$STAGE"/$(basename {output.noncpgfile:q}) \
            --bed-methyl "$STAGE"/{wildcards.sample:q} \
            --tabix \
            --threads {threads} \
            {params.flags:q} \
            {input.bcf:q}

        # Each .md5 names its file as mextr was given it.
        OUT=$(dirname {output.cpgfile:q})
        sed -i "s|$STAGE/|$OUT/|" "$STAGE"/*.md5
        mv "$STAGE"/* "$OUT"/
        rmdir "$STAGE"
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
    params:
        tmpdir=config["tmpdir"],
    log:
        f"{LOGDIR}/extract_signal/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/extract_signal/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        mkdir -p {params.tmpdir:q}
        python {workflow.basedir:q}/scripts/extract_signal.py \
            --bed-gz {input.bed_gz:q} \
            --chromsizes {input.chromsizes:q} \
            --output-dir $(dirname {output.bigwigs[0]:q}) \
            --barcode {wildcards.sample:q} \
            --tmpdir {params.tmpdir:q}
        """


rule get_methylc:
    input:
        chh_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_chh.bed.gz",
        chg_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_chg.bed.gz",
        cpg_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_cpg.bed.gz",
    output:
        methylc_bed_gz=f"{OUTDIR}/{{sample}}/{{sample}}_methylc.bed.gz",
        methylc_tbi=f"{OUTDIR}/{{sample}}/{{sample}}_methylc.bed.gz.tbi",
    params:
        tmpdir=config["tmpdir"],
    log:
        f"{LOGDIR}/get_methylc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/get_methylc/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        mkdir -p {params.tmpdir:q}
        python {workflow.basedir:q}/scripts/get_methylc.py \
            --barcode {wildcards.sample:q} \
            --chh {input.chh_bed_gz:q} \
            --chg {input.chg_bed_gz:q} \
            --cpg {input.cpg_bed_gz:q} \
            --threads {threads} \
            --outdir $(dirname {output.methylc_bed_gz:q}) \
            --tmpdir {params.tmpdir:q}
        """
