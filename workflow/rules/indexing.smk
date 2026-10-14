rule indexing:
    input:
        fasta=REFERENCE,
    output:
        fasta=f"{REFDIR}/reference.fasta",
        fasta_fai=f"{REFDIR}/reference.fasta.fai",
        chromsizes=f"{REFDIR}/chromsizes.tsv",
        contig_md5=f"{REFDIR}/reference.contig_md5",
    params:
        dict_flags=lambda w: dict_flags(),
    log:
        f"{LOGDIR}/indexing.txt",
    benchmark:
        f"{BENCHDIR}/indexing.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        if [[ {input.fasta:q} == *.gz ]]; then
            unpigz -p {threads} -c {input.fasta:q} > {output.fasta:q}
        else
            cp {input.fasta:q} {output.fasta:q}
        fi
        samtools faidx {output.fasta:q}
        cut -f1,2 {output.fasta_fai:q} > {output.chromsizes:q}

        # gemBS's contig_md5 format: the @SQ fields without the record type and SN:.
        samtools dict {params.dict_flags:q} {output.fasta:q} \
            | sed -n 's/^@SQ\tSN://p' > {output.contig_md5:q}

        sync {output:q}
        """


rule gem_indexer:
    input:
        fasta=f"{REFDIR}/reference.fasta",
    output:
        gem_index=f"{REFDIR}/reference.gem",
        gem_info=f"{REFDIR}/reference.info",
    params:
        sampling_rate=lambda w: config["index"]["sampling_rate"],
        tmpdir=config["tmpdir"],
    log:
        f"{LOGDIR}/gem_indexer.txt",
    benchmark:
        f"{BENCHDIR}/gem_indexer.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # gem-indexer appends file names to --tmp-folder without a slash.
        mkdir -p {params.tmpdir:q}
        gem-indexer \
            -i {input.fasta:q} \
            -o $(dirname {output.gem_index:q})/reference \
            --bisulfite-index \
            --text-sampling-rate {params.sampling_rate:q} \
            --threads {threads} \
            --tmp-folder {params.tmpdir:q}/

        sync {output:q}
        """


rule bwamem3_index:
    input:
        fasta=f"{REFDIR}/reference.fasta",
    output:
        index=multiext(
            f"{REFDIR}/reference.fasta",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
            ".meth.amb",
            ".meth.ann",
            ".meth.bwt.2bit.64",
            ".meth.pac",
        ),
        # The converted FASTA of the seed index, which mem --meth does not read.
        meth_fasta=temp(f"{REFDIR}/reference.fasta.meth.fa"),
    log:
        f"{LOGDIR}/bwamem3_index.txt",
    benchmark:
        f"{BENCHDIR}/bwamem3_index.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # Scratch beside the index, off the node's shared /tmp.
        bwa-mem3 index \
            --meth \
            -t {threads} \
            --tmp-dir $(dirname {input.fasta:q}) \
            {input.fasta:q}

        sync {output:q}
        """
