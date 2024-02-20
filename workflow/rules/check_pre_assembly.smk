rule check_foreign_contamination:
    input:
        "results/kraken/{sample}.bracken",
    output:
        temp("results/checks/{sample}/foreign_contamination.tsv"),
    params:
        fraction_threshold=config["foreign_contamination"]["abundance_check_fraction"],
    log:
        "logs/checks/foreign_contamination/{sample}.log",
    localrule: True
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/genera_check.py"


rule check_number_of_bases:
    input:
        "results/reads/trimmed/{sample}.qc.txt",
    output:
        temp("results/checks/{sample}/basepairs_for_assembly.tsv"),
    params:
        min_bp=config["min_basepairs_for_assembly"],
    log:
        "logs/checks/check_number_of_bases/{sample}.log",
    localrule: True
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/read_bp_check.py"


checkpoint pre_assembly_QC:
    input:
        bracken="results/checks/{sample}/foreign_contamination.tsv",
        bps="results/checks/{sample}/basepairs_for_assembly.tsv",
    output:
        "results/checks/{sample}/pre_assembly_summary.tsv",
    log:
        "logs/checks/pre_assembly/{sample}.log",
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    shell:
        "cat {input} > {output} 2>&1"
