rule amrfinder__download_db:
    output:
        db=directory(config["amrfinder_db_dir"]),
    params:
        db_parent=lambda wildcards, output: os.path.dirname(output.db),
    conda:
        "../envs/amrfinder.yaml"
    log:
        os.path.join(os.path.dirname(config["amrfinder_db_dir"]), "logs", "download"),
    shell:
        "(mkdir -p {params.db_parent} && amrfinder_update -d {params.db_parent}) > {log} 2>&1"


rule amrfinder__call:
    input:
        contigs=infer_assembly_fasta,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
        db=config["amrfinder_db_dir"],
    output:
        tsv="results/amr_detect/{sample}/amrfinder.tsv",
    params:
        organism_arg=get_organism_for_amrfinder,
    threads: min(config["threads"]["amrfinder"], config["max_threads"])
    conda:
        "../envs/amrfinder.yaml"
    log:
        "logs/amr_detect/amrfinder/{sample}.log",
    shell:
        "amrfinder -d {input.db} --nucleotide {input.contigs} --threads {threads} {params.organism_arg} -o {output.tsv} > {log} 2>&1"


rule mlst__call:
    input:
        contigs=infer_assembly_fasta,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
    output:
        "results/amr_detect/{sample}/mlst.tsv",
    params:
        scheme_arg=get_taxonomy_for_mlst,
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/amr_detect/mlst/{sample}.log",
    shell:
        "mlst {input.contigs} {params.scheme_arg} > {output} 2> {log}"


rule abricate__call:
    input:
        infer_assembly_fasta,
    output:
        "results/amr_detect/{sample}/abricate.tsv",
    params:
        abricate_db=config["abricate"]["db"],
        min_identity=config["abricate"]["min_identity"],
        min_coverage=config["abricate"]["min_coverage"],
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/amr_detect/abricate/{sample}.log",
    shell:
        "abricate --db {params.abricate_db} --minid {params.min_identity} --mincov {params.min_coverage}"
        " {input} > {output} 2> {log}"


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


rule spatyper__database_download:
    output:
        repeats=protected(os.path.join(config["spatyper_db_dir"], "sparepeats.fasta")),
        order=protected(os.path.join(config["spatyper_db_dir"], "spatypes.txt")),
    params:
        db_dir=lambda wildcards, output: os.path.dirname(output.repeats),
    localrule: True
    conda:
        "../envs/spatyper.yaml"
    log:
        os.path.join(os.path.join(config["spatyper_db_dir"], "logs", "download.log")),
    script:
        "../scripts/spatyper_db_download.py"


rule spatyper__call:
    input:
        fasta=infer_assembly_fasta,
        repeats=os.path.join(config["spatyper_db_dir"], "sparepeats.fasta"),
        order=os.path.join(config["spatyper_db_dir"], "spatypes.txt"),
    output:
        "results/amr_detect/{sample}/spa_typer.tsv",
    conda:
        "../envs/spatyper.yaml"
    log:
        "logs/amr_detect/spa_typer/{sample}.log",
    shell:
        "spaTyper -f {input.fasta} -r {input.repeats} -o {input.order} --output {output} > {log} 2>&1"


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


rule sccmec__download_db:
    output:
        db_proteins=protected(os.path.join(config["SCCmec_db_dir"], "proteins.fasta")),
        db_subtypes=protected(os.path.join(config["SCCmec_db_dir"], "subtypes.fasta")),
        db_primers=protected(os.path.join(config["SCCmec_db_dir"], "primers.fasta")),
    params:
        repo="https://github.com/staphopia/staphopia-sccmec/archive/refs/heads/master.zip",
        repo_db_name="staphopia-sccmec-master/share/staphopia-sccmec/data",
        db_dir=lambda wildcards, output: os.path.dirname(output.db_proteins),
    conda:
        "../envs/curl_with_unzip.yaml"
    log:
        os.path.join(config["SCCmec_db_dir"], "logs", "download.log"),
    script:
        "../scripts/sccmec_download_db.py"


rule sccmec__call:
    input:
        assembly=infer_assembly_fasta,
        db_proteins=os.path.join(config["SCCmec_db_dir"], "proteins.fasta"),
        db_subtypes=os.path.join(config["SCCmec_db_dir"], "subtypes.fasta"),
        db_primers=os.path.join(config["SCCmec_db_dir"], "primers.fasta"),
    output:
        "results/amr_detect/{sample}/SCCmec.tsv",
    params:
        ext=lambda wildcards, input: os.path.splitext(input[0])[1],
        db_dir=lambda wildcards, input: os.path.dirname(input.db_proteins),
    conda:
        "../envs/sccmec.yaml"
    log:
        "logs/amr_detect/{sample}/sccmec.log",
    shell:
        "staphopia-sccmec --assembly {input.assembly} --sccmec {params.db_dir} --ext {params.ext} > {output} 2>{log}"


rule mob_suite__download_db:
    output:
        db=protected(os.path.join(config["ncbi_plasmid_db_dir"], "ncbi_plasmid_full_seqs.fas.msh")),
    params:
        db_dir=lambda wildcards, output: os.path.dirname(output.db),
    conda:
        "../envs/mob_suite.yaml"
    log:
        "logs/plasmids/download_db.log",
    shell:
        "mob_init --database_directory {params.db_dir} > {log} 2>&1"


rule mob_suite__typer:
    input:
        fasta=infer_assembly_fasta,
        db=os.path.join(config["ncbi_plasmid_db_dir"], "ncbi_plasmid_full_seqs.fas.msh"),
    output:
        "results/plasmids/{sample}/mob_typer.txt",
    params:
        db_dir=lambda wildcards, input: os.path.dirname(input.db),
    conda:
        "../envs/mob_suite.yaml"
    log:
        "logs/plasmids/{sample}/mob_typer.log",
    shell:
        "mob_typer --infile {input.fasta} --database_directory {params.db_dir} --out_file {output} > {log} 2>&1"


rule sistr_cmd__call:
    input:
        infer_assembly_fasta,
    output:
        serovar="results/amr_detect/{sample}/sistr_serovar.csv",
        alleles="results/amr_detect/{sample}/sistr/allele-results.json",
        cgmlst="results/amr_detect/{sample}/sistr/cgmlst-profiles.csv",
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.alleles),
    conda:
        "../envs/sistr_cmd.yaml"
    log:
        "logs/amr_detect/sistr/{sample}.log",
    shell:
        "(mkdir -p {params.out_dir} && sistr --qc --output-format tab --output-prediction {output.serovar}"
        " --alleles-output {output.alleles} --cgmlst-profiles {output.cgmlst}) > {log} 2>&1"


rule rgi_download_db:
    output:
        json=protected(os.path.join(config["rgi_db_dir"], "card.json")),
    params:
        db_dir=lambda wildcards, output: os.path.dirname(output.json),
        db_url="https://card.mcmaster.ca/latest/data",
    conda:
        "../envs/curl.yaml"
    log:
        os.path.join(config["rgi_db_dir"], "logs", "download.log"),
    shell:
        "mkdir -p {params.db_dir} && (curl -SL {params.db_url} | tar xjvf - -C {params.db_dir}) > {log} 2>&1"


rule rgi_load_db:
    input:
        json=os.path.join(config["rgi_db_dir"], "card.json"),
    output:
        loaded_db=temp(directory("localDB")),
    conda:
        "../envs/rgi.yaml"
    log:
        "logs/amr_detect/rgi_db_load.log",
    shell:
        "(rm -rf {output.loaded_db} && rgi load --card_json {input.json} --local) > {log} 2>&1"


rule rgi_call:
    input:
        json=os.path.join(config["rgi_db_dir"], "card.json"),
        loaded_db="localDB/",
        assembly=infer_assembly_fasta,
    output:
        txt="results/amr_detect/{sample}/rgi_main.txt",
        json="results/amr_detect/{sample}/rgi_main.json",
    params:
        extra="--include_nudge",
        out_prefix=lambda wildcards, output: os.path.splitext(output.txt)[0],
    conda:
        "../envs/rgi.yaml"
    threads: min(config["threads"]["amrfinder"], config["max_threads"])
    log:
        "logs/amr_detect/rgi/{sample}.log",
    shell:
        "rgi main --input_sequence {input.assembly} --local --clean {params.extra}"
        " --output_file {params.out_prefix} --num_threads {threads} --split_prodigal_jobs > {log} 2>&1"
