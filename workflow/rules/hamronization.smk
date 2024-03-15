rule abricate__version_tool:
    output:
        "results/.versions/abricate_tool.txt",
    localrule: True
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/versions/abricate_tool.log",
    shell:
        "(abricate --version | tr -s '[:blank:]' '\n' | tail -1 > {output}) 2> {log}"


rule amrfinder__version_tool:
    output:
        "results/.versions/amrfinder_tool.txt",
    localrule: True
    conda:
        "../envs/amrfinder.yaml"
    log:
        "logs/versions/amrfinder_tool.log",
    shell:
        "amrfinder --version > {output} 2> {log}"


rule rgi__version_tool:
    output:
        "results/.versions/rgi_tool.txt",
    localrule: True
    conda:
        "../envs/rgi.yaml"
    log:
        "logs/versions/rgi_tool.log",
    shell:
        "rgi main --version > {output} 2> {log}"


rule abricate__version_db:
    input:
        db=os.path.join(config["abricate"]["db_dir"], "sequences"),
    output:
        "results/.versions/abricate_db.txt",
    params:
        db_name=lambda wildcards, input: os.path.basename(os.path.dirname(input.db)),
        db_dir=lambda wildcards, input: os.path.dirname(os.path.dirname(input.db)),
    localrule: True
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/versions/abricate_db.log",
    shell:
        "(abricate --list --datadir {params.db_dir} | grep {params.db_name}"
        " | tr -s '[:blank:]' '\n' | tail -1 > {output}) 2> {log}"


rule amrfinder__version_db:
    input:
        db=config["amrfinder_db_dir"],
    output:
        "results/.versions/amrfinder_db.txt",
    conda:
        "../envs/amrfinder.yaml"
    localrule: True
    log:
        "logs/versions/amrfinder_db.log",
    shell:
        "(amrfinder -V -d {input.db} | grep 'Database version:'"
        " | tr ':' '\n' | tail -1 | tr -d ' ' > {output}) 2> {log}"


rule rgi__version_db:
    input:
        db="localDB/",
    output:
        "results/.versions/rgi_db.txt",
    params:
        json=lambda wildcards, input: os.path.join(input.db, "loaded_databases.json"),
    localrule: True
    conda:
        "../envs/python.yaml"
    log:
        "logs/versions/rgi_db.log",
    script:
        "../scripts/rgi_get_version_db.py"


rule hamronize__rgi:
    input:
        result="results/amr_detect/{sample}/rgi_main.txt",
        version="results/.versions/rgi_db.txt",
        db_version="results/.versions/rgi_tool.txt",
    output:
        tsv="results/hamronization/rgi/{sample}.tsv",
    localrule: True
    params:
        sample_name=lambda wildcards: wildcards.sample,
        tool="rgi",
    log:
        "logs/hamronization/rgi/{sample}.log",
    conda:
        "../envs/hamronize.yaml"
    script:
        "../scripts/hamronize.py"


rule hamronize__amrfinder:
    input:
        result="results/amr_detect/{sample}/amrfinder.tsv",
        version="results/.versions/amrfinder_db.txt",
        db_version="results/.versions/amrfinder_tool.txt",
    output:
        tsv="results/hamronization/amrfinder/{sample}.tsv",
    localrule: True
    params:
        sample_name=lambda wildcards: wildcards.sample,
        tool="amrfinder",
    log:
        "logs/hamronization/amrfinder/{sample}.log",
    conda:
        "../envs/hamronize.yaml"
    script:
        "../scripts/hamronize.py"


rule hamronize__abricate:
    input:
        result="results/amr_detect/{sample}/abricate.tsv",
        version="results/.versions/abricate_db.txt",
        db_version="results/.versions/abricate_tool.txt",
    output:
        tsv="results/hamronization/abricate/{sample}.tsv",
    localrule: True
    params:
        tool="abricate",
    log:
        "logs/hamronization/abricate/{sample}.log",
    conda:
        "../envs/hamronize.yaml"
    script:
        "../scripts/hamronize.py"


rule hamronize__summarize:
    input:
        infer_amr_detection_results_for_harmonize,
    output:
        tsv="results/hamronization/summary.{ext}",
    params:
        type_format=lambda wildcards, output: "tsv" if wildcards.ext == "tsv" else "interactive",
    log:
        "logs/hamronization/summary_{ext}.log",
    conda:
        "../envs/hamronize.yaml"
    shell:
        "hamronize summarize -o {output} -t {params.type_format} {input} > {log} 2>&1"


rule request_hamronize_summary:
    input:
        request_hamronize_or_nothing,
    output:
        "results/hamronization/hamronization_requested.txt",
    params:
        value=lambda wildcards, input: "REQUESTED" if len(input) > 0 else "NOT_REQUESTED",
    conda:
        "../envs/coreutils.yaml"
    log:
        "logs/checks/request_hamronize_summary.log",
    shell:
        "echo {params.value} > {output} 2> {log}"
