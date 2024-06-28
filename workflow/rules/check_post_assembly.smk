checkpoint checkpoint_assembly_construction:
    input:
        "results/assembly/{sample}/assembly.gfa",
    output:
        "results/checks/{sample}/assembly_constructed.tsv",
    params:
        happy_msg="\t".join(["PASS", "assembly_construction", "success", "Assembly is not empty"]),
        sad_msg="\t".join(["FAIL", "assembly_construction", "failure", "Assembly construction failed"]),
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    log:
        "logs/checks/assembly_constructed/{sample}.log",
    shell:
        "(([ -s {input} ] && echo {params.happy_msg:q}) || echo {params.sad_msg:q}) > {output} 2> {log}"


rule log_that_assembly_was_not_requested:
    input:
        infer_assembly_fasta,
    output:
        "results/checks/{sample}/check_skipping.tsv",
    conda:
        "../envs/coreutils.yaml"
    params:
        message="\t".join(["PASS", "assembly_not_requested", "true", "Assembly provided as input"]),
    log:
        "logs/checks/assembly/{sample}.log",
    localrule: True
    shell:
        "echo {params.message:q} > {output} 2> {log}"


rule check_assembly_quality:
    input:
        tsv="results/assembly/{sample}/bandage/bandage.info",
        svg="results/assembly/{sample}/bandage/bandage.svg",
        stats="results/assembly/{sample}/seqkit_stats.tsv",
    output:
        temp("results/checks/{sample}/assembly_quality.tsv"),
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
        vcf="results/self_contamination/{sample}/filtered.vcf",
        txt="results/self_contamination/{sample}/markdup_isize.txt", # NOTE just to request
        stats_orig="results/self_contamination/{sample}/mapped_stats.txt", # NOTE just to request
        stats_dup="results/self_contamination/{sample}/markdup_stats.txt", # NOTE just to request
    output:
        temp("results/checks/{sample}/self_contamination_check.tsv"),
    params:
        max_rows=config["self_contamination"]["max_ambiguous_rows"],
        check_level=config["self_contamination"]["check_level"],
    log:
        "logs/checks/self_contamination/{sample}.log",
    conda:
        "../envs/python.yaml"
    localrule: True
    script:
        "../scripts/self_contamination.py"


rule check_coverage_from_qualimap:
    input:
        "results/self_contamination/{sample}/markdup/bamqc",
    output:
        temp("results/checks/{sample}/coverage_check.tsv"),
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


checkpoint checkpoint_request_post_assembly_checks_if_relevant:
    input:
        infer_relevant_checks,
    output:
        "results/checks/{sample}/qc_summary.tsv",
    log:
        "logs/checks/summary/{sample}.log",
    params:
        header="result\tparameter\tvalue\tcomment",
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    shell:
        "(echo -e {params.header:q} && cat {input}) > {output} 2> {log}"


rule request_all_relevant_outputs:
    input:
        infer_outputs_for_sample_as_list,
    output:
        temp("results/checks/{sample}/.final_results_requested.tsv"),
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    log:
        "logs/checks/assembly/{sample}.log",
    shell:
        "echo {input} > {output} 2> {log}"
