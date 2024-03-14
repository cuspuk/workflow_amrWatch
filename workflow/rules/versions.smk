rule abricate__version_tool:
    output:
        "results/.versions/abricate_tool.txt",
    localrule: True
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/versions/abricate_tool.log",
    shell:
        "abricate --version > {output} 2> {log}"


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
    """
    Produces a file with the version of the abricate database.
    Example:
    vfdb  4366  nucl  2024-Mar-14
    """
    input:
        db=protected(os.path.join(config["abricate"]["db_dir"], "sequences")),
    output:
        "results/.versions/abricate_db.tsv",
    params:
        db_name=lambda wildcards, input: os.path.basename(os.path.dirname(input.db)),
        db_dir=lambda wildcards, input: os.path.dirname(os.path.dirname(input.db)),
    localrule: True
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/versions/abricate_db.log",
    shell:
        "abricate --list --datadir {params.db_dir} | grep {params.db_name}"
        " | | tr -s '[:blank:]' '\n' | tail -1 > {output} 2> {log}"


rule amrfinder__version_db:
    input:
        db=directory(config["amrfinder_db_dir"]),
    output:
        "results/.versions/amrfinder_db.txt",
    conda:
        "../envs/amrfinder.yaml"
    log:
        "logs/versions/amrfinder_db.log",
    shell:
        "(amrfinder -V -d {input.db} | grep 'Database version:'"
        " | tr ':' '\n' | tail -1 | tr -d ' ') > {output} 2> {log}"


rule rgi__version_db:
    input:
        db="localDB/",
    output:
        "results/.versions/rgi_db.txt",
    params:
        json=lambda wildcards, input: os.path.join(input.db, "loaded_databases.json"),
    conda:
        "../envs/python.yaml"
    log:
        "logs/versions/rgi_db.log",
    script:
        "../scripts/rgi__version_db.py"
