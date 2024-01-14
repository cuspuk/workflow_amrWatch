rule abritamr__call:
    input:
        "results/assembly/{sample}/assembly.fasta",
    output:
        tsv="results/amr_detect/{sample}/amrfinder.tsv",
        dir=directory("results/amr_detect/{sample}/amrfinder"),
    params:
        identity=0.9,  # TODO
        tmp_output=lambda wildcards, output: os.path.join(output.dir, "amrfinder.out"),
    threads: min(config["threads"]["abritamr"], config["max_threads"])
    conda:
        "../envs/abritamr.yaml"
    log:
        "logs/amr_detect/abritamr/{sample}.log",
    shell:
        "(abritamr run --contigs {input} --prefix {output.dir} --jobs {threads} --identity {params.identity}"
        " && mv {params.tmp_output} {output.tsv} ) > {log} 2>&1"


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
        "mlst {input} > {output} 2> {log}"


rule abricate__call:
    input:
        "results/assembly/{sample}/assembly.fasta",
    output:
        "results/amr_detect/{sample}/abricate.tsv",
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/amr_detect/mlst/{sample}.log",
    shell:
        "abricate {input} > {output} 2> {log}"
