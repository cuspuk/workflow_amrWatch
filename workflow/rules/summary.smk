
rule summary__results_per_sample:
    input:
        unpack(infer_results_to_summarize_for_sample),
    output:
        tsv="results/summary/per_sample/{sample}.tsv",
    params:
        delimiter="\t",
        sample_name=lambda wildcards: wildcards.sample,
        amrfinder_uniq_tag="__",
    localrule: True
    conda:
        "../envs/python.yaml"
    log:
        "logs/summary/summary_results/{sample}.log",
    script:
        "../scripts/summary.py"


rule merge__summary_results_per_sample:
    input:
        tsvs=expand("results/summary/per_sample/{sample}.tsv", sample=get_sample_names()),
    output:
        tsv="results/summary/summary.tsv",
        amrfinder_file="results/summary/amrfinder_output.tsv",
    params:
        delimiter="\t",
        amrfinder_uniq_tag="__",
        nan_value="NaN",
    localrule: True
    conda:
        "../envs/python.yaml"
    log:
        "logs/summary/merge_summaries.log",
    script:
        "../scripts/merge_summaries.py"
