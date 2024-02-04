checkpoint assembly_constructed:
    input:
        "results/assembly/{sample}/assembly.gfa",
    output:
        "results/checks/{sample}/assembly_constructed.txt",
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    log:
        "logs/checks/assembly_constructed/{sample}.log",
    shell:
        "(([ -s {input} ] && echo 'PASS: Assembly is not empty') || echo 'FAIL: Assembly construction failed') > {output} 2> {log}"


rule assembly_not_requested:
    output:
        "results/checks/{sample}/check_skipping.txt",
    conda:
        "../envs/coreutils.yaml"
    log:
        "logs/checks/assembly/{sample}.log",
    localrule: True
    shell:
        "echo 'PASS: Assembly provided as input' > {output} 2> {log}"


rule check_assembly_quality:
    input:
        tsv="results/assembly/{sample}/bandage/bandage.info",
        svg="results/assembly/{sample}/bandage/bandage.svg",
    output:
        "results/checks/{sample}/assembly_quality.txt",
    params:
        max_dead_ends=config["assembly__unicycler"]["max_dead_ends"],
        max_contigs=config["assembly__unicycler"]["max_contigs"],
        min_length_in_bp=config["assembly__unicycler"]["min_length_in_bp"],
        max_length_in_bp=config["assembly__unicycler"]["max_length_in_bp"],
    conda:
        "../envs/python.yaml"
    localrule: True
    log:
        "logs/checks/assembly_quality/{sample}.log",
    script:
        "../scripts/check_assembly_quality.py"


rule check_self_contamination:
    input:
        "results/self_contamination/{sample}/filtered.vcf",
    output:
        "results/checks/{sample}/self_contamination_check.txt",
    params:
        max_rows=config["self_contamination"]["max_ambiguous_rows"],
    log:
        "logs/checks/self_contamination/{sample}.log",
    conda:
        "../envs/grep.yaml"
    localrule: True
    script:
        "../scripts/self_contamination.sh"


rule check_coverage_from_qualimap:
    input:
        "results/self_contamination/{sample}/markdup/bamqc",
    output:
        "results/checks/{sample}/coverage_check.txt",
    params:
        genome_results_file=lambda wildcards, input: os.path.join(input[0], "genome_results.txt"),
        warn_threshold=config["coverage_check"]["warn_threshold"],
        fail_threshold=config["coverage_check"]["fail_threshold"],
    log:
        "logs/checks/coverage_check/{sample}.log",
    conda:
        "../envs/python.yaml"
    localrule: True
    script:
        "../scripts/coverage_check.py"


rule check_foreign_contamination:
    input:
        "results/kraken/{sample}.bracken",
    output:
        "results/checks/{sample}/foreign_contamination.txt",
    params:
        fraction_threshold=config["foreign_contamination"]["abundance_check_fraction"],
    log:
        "logs/checks/foreign_contamination/{sample}.log",
    localrule: True
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/genera_check.py"


checkpoint summary_all_checks:
    input:
        infer_relevant_checks,
    output:
        "results/checks/{sample}/summary.txt",
    log:
        "logs/checks/summary/{sample}.log",
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    shell:
        "cat {input} > {output} 2>&1"


rule get_final_outputs:
    input:
        infer_outputs_for_sample,
    output:
        temp("results/checks/{sample}/.final_results_requested.txt"),
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    log:
        "logs/checks/assembly/{sample}.log",
    shell:
        "echo {input} > {output} 2> {log}"
