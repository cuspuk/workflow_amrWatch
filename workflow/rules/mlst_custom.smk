rule mlst__call__custom:
    input:
        contigs=infer_assembly_fasta,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
        blast_mlst=os.path.join(config["mlst_custom"]["db_dir"], "blast", "mlst.fa"),
        pubmlst=os.path.join(config["mlst_custom"]["db_dir"], "pubmlst"),
    output:
        "results/amr_detect/{sample}/mlst_custom.tsv"
    params:
        scheme_arg=get_custom_scheme_for_mlst,
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/amr_detect/mlst_custom/{sample}.log",
    shell:
        "mlst {input.contigs} --blastdb {input.blast_mlst} --datadir {input.pubmlst} {params.scheme_arg} --legacy > {output} 2> {log}"


rule mlst_make_custom_blast_db:
    input:
        mlst_dir=os.path.join(config["mlst_custom"]["db_dir"], "pubmlst"),
        schemas=infer_custom_pubmlst_schemas,
    output:
        blast_dir=directory(os.path.join(config["mlst_custom"]["db_dir"], "blast")),
        blast_file=os.path.join(config["mlst_custom"]["db_dir"], "blast", "mlst.fa"),
    localrule: True
    conda:
        "../envs/mlst.yaml"
    log:
        os.path.join(config["mlst_custom"]["db_dir"], "logs", "make_blast_db.log"),
    shell:
        """
        (mkdir -p {output.blast_dir}
        rm -f {output.blast_file}

        for N in $(find {input.mlst_dir} -mindepth 1 -maxdepth 1 -type d); do
        SCHEME=$(basename $N)
        echo "Adding: $SCHEME"
        cat "{input.mlst_dir}"/$SCHEME/*.tfa \
            | grep -v 'not a locus'  \
            | sed -e "s/^>/>$SCHEME./" \
            >> {output.blast_file}
        done

        makeblastdb -hash_index -in {output.blast_file} -dbtype nucl -title "PubMLST" -parse_seqids) > {log} 2>&1
        """