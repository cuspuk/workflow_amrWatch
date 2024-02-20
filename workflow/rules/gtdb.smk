rule gtdbtk__classify:
    input:
        assembly=infer_assembly_fasta,
        gtdb=os.path.join(config["gtdb_dirpath"], "db"),
    output:
        gtdb_tsv="results/taxonomy/{sample}/classify/gtdbtk.bac120.summary.tsv",
    params:
        assembly_dir=lambda wildcards, input: os.path.dirname(input.assembly),
        out_dir=lambda wildcards, output: os.path.dirname(os.path.dirname(output.gtdb_tsv)),
    threads: min(config["threads"]["gtdb__classify"], config["max_threads"])
    conda:
        "../envs/gtdbtk.yaml"
    resources:
        mem_mb=get_mem_mb_for_gtdb,
    log:
        "logs/taxonomy/gtdb_classify/{sample}.log",
    shell:
        "(rm -rf {params.out_dir} && export GTDBTK_DATA_PATH={input.gtdb:q} "
        " && echo '{input.assembly}\t{wildcards.sample}' > $TMPDIR/batchfile_{wildcards.sample}.txt"
        " && gtdbtk classify_wf --batchfile $TMPDIR/batchfile_{wildcards.sample}.txt"
        " --cpus {threads} --tmpdir $TMPDIR --extension fasta --out_dir {params.out_dir} --mash_db {input.gtdb}) > {log} 2>&1"


rule gtdbtk__download_metadata:
    output:
        protected("{gtdb_dir}/bac120_metadata.tsv"),
    params:
        gzipped=lambda wildcards, output: f"{output}.gz",
    conda:
        "../envs/coreutils.yaml"
    log:
        "{gtdb_dir}/download_metadata.log",
    shell:
        "(wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz -O {params.gzipped} && gzip -d {params.gzipped}) > {log} 2>&1"


rule gtdbtk__convert_to_ncbi:
    input:
        metadata=os.path.join(config["gtdb_dirpath"], "bac120_metadata.tsv"),
        gtdb_tsv="results/taxonomy/{sample}/classify/gtdbtk.bac120.summary.tsv",
    output:
        "results/taxonomy/{sample}/ncbi_taxa.tsv",
    params:
        gtdb_parent_dir=lambda wildcards, input: os.path.dirname(os.path.dirname(input.gtdb_tsv)),
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdb_convert_to_ncbi/{sample}.log",
    script:
        "../scripts/gtdbtk_ncbi_convert.py"


checkpoint checkpoint_parse_taxa_gtdbtk:
    input:
        gtdb_tsv="results/taxonomy/{sample}/classify/gtdbtk.bac120.summary.tsv",
        ncbi_taxa="results/taxonomy/{sample}/ncbi_taxa.tsv",  # not a required dependency. Specified to simplify gathering of results.
    output:
        "results/taxonomy/{sample}/parsed_taxa.txt",
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    log:
        "logs/taxonomy/parse_taxa/{sample}.log",
    shell:
        '(cut -f2 {input.gtdb_tsv} | tail -n 1 | sed -e "s/.*;s__//") > {output} 2> {log}'
