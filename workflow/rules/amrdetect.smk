rule abritamr__call:
    input:
        "results/assembly/{sample}/assembly.fasta",
    output:
        "results/amr_detect/{sample}/amrfinder",
    params:
        identity=0.9,  # TODO
    threads: min(config["threads"]["abritamr"], config["max_threads"])
    conda:
        "../envs/abritamr.yaml"
    log:
        "logs/amr_detect/abritamr/{sample}.log",
    shell:
        "abritamr run --contigs {input} --prefix {output} --jobs {threads} --identity {params.identity} > {log} 2>&1"


rule mlst__call:
    input:
        "results/assembly/{sample}/assembly.fasta",
    output:
        "results/amr_detect/{sample}/mlst.tsv",
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/amr_detect/mlst/{sample}.log",
    shell:
        "mlst {input} > {output} 2>&1"
