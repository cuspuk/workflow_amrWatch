rule amrfinder__call:
    output:
        db=directory(config["amrfinder_db_dir"]),
    conda:
        "../envs/amrfinder.yaml"
    log:
        os.path.join(config["amrfinder_db_dir"], "logs", "download"),
    shell:
        "amrfinder_update -d {output.db} > {log} 2>&1"


rule amrfinder__call:
    input:
        contigs=infer_assembly_fasta,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
        db=config["amrfinder_db_dir"],
    output:
        tsv="results/amr_detect/{sample}/amrfinder.tsv",
    params:
        organism=get_organism_for_amrfinder,
    threads: min(config["threads"]["amrfinder"], config["max_threads"])
    conda:
        "../envs/amrfinder.yaml"
    log:
        "logs/amr_detect/amrfinder/{sample}.log",
    shell:
        "amrfinder -d {input.db} --nucleotide {input.contigs} --threads {threads} --organism {params.organism} -o {output.tsv} > {log} 2>&1"


rule mlst__call:
    input:
        contigs=infer_assembly_fasta,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
    output:
        "results/amr_detect/{sample}/mlst.tsv",
    params:
        scheme=get_taxonomy_for_mlst,
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/amr_detect/mlst/{sample}.log",
    shell:
        "mlst {input.contigs} --scheme {params.scheme} > {output} 2> {log}"


rule abricate__call:
    input:
        infer_assembly_fasta,
    output:
        "results/amr_detect/{sample}/abricate.tsv",
    params:
        abricate_db=config["abricate_db"],
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/amr_detect/abricate/{sample}.log",
    shell:
        "abricate --db {params.abricate_db} {input} > {output} 2> {log}"


rule kleborate__call:
    input:
        infer_assembly_fasta,
    output:
        "results/amr_detect/{sample}/kleborate.tsv",
    conda:
        "../envs/kleborate.yaml"
    log:
        "logs/amr_detect/kleborate/{sample}.log",
    shell:
        "kleborate --all -o {output} -a {input} > {log} 2>&1"


rule spatyper__call:
    input:
        infer_assembly_fasta,
    output:
        "results/amr_detect/{sample}/spa_typer.tsv",
    conda:
        "../envs/spatyper.yaml"
    log:
        "logs/amr_detect/spa_typer/{sample}.log",
    shell:
        "spaTyper -f {input} --output {output} > {log} 2>&1"


rule etoki__call:
    input:
        infer_assembly_fasta,
    output:
        "results/amr_detect/{sample}/etoki_ebeis.tsv",
    conda:
        "../envs/etoki.yaml"
    log:
        "logs/amr_detect/etoki/{sample}.log",
    shell:
        "EToKi.py EBEis -q {input} > {output} 2>{log}"
