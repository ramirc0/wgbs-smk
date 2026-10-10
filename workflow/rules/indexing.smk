rule indexing:
    input:
        fasta=REFERENCE,
    output:
        fasta=f"{REFDIR}/reference.fasta",
        fasta_fai=f"{REFDIR}/reference.fasta.fai",
        chromsizes=f"{REFDIR}/chromsizes.tsv",
        sam_dict=f"{REFDIR}/sam_dict.tsv",
    params:
        stage=lambda w: f"{config['tmpdir']}/{REFDIR}",
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

        OUTPUT_DIR={params.stage:q}

        mkdir -p "$OUTPUT_DIR"
        cp {input.fasta:q} "$OUTPUT_DIR"
        FASTA=$OUTPUT_DIR/$(basename {input.fasta:q})
        if [[ "$FASTA" == *.gz ]]; then
            unpigz -c "$FASTA" > $OUTPUT_DIR/$(basename {output.fasta:q})
        else
            mv "$FASTA" $OUTPUT_DIR/$(basename {output.fasta:q})
        fi
        rm -rf "$FASTA"

        samtools faidx -@ {threads} $OUTPUT_DIR/$(basename {output.fasta:q})
        cut -f1,2 $OUTPUT_DIR/$(basename {output.fasta:q}).fai > $OUTPUT_DIR/$(basename {output.chromsizes:q})
        samtools dict \
            -o $OUTPUT_DIR/$(basename {output.sam_dict:q}) \
            {params.dict_flags:q} \
            $OUTPUT_DIR/$(basename {output.fasta:q})

        cp -t $(dirname {output.fasta:q}) "$OUTPUT_DIR"/*
        rm -rf "$OUTPUT_DIR"

        echo "Done!"
        """


rule gem_indexer:
    input:
        fasta=f"{REFDIR}/reference.fasta",
    output:
        gem_index=f"{REFDIR}/reference.gem",
    params:
        sampling_rate=lambda w: config["index"]["sampling_rate"],
    log:
        f"{LOGDIR}/gem_indexer.txt",
    benchmark:
        f"{BENCHDIR}/gem_indexer.tsv",
    shadow:
        "copy-minimal"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        gem-indexer \
            -i {input.fasta:q} \
            -o $(dirname {output.gem_index:q})/reference \
            --bisulfite-index \
            --text-sampling-rate {params.sampling_rate:q} \
            --threads {threads}

        echo "Done!"
        """
