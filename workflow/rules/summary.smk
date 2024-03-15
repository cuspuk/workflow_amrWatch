
rule summary__results_per_sample:
    input:
        infer_results_to_summarize_for_sample,
    output:
        tsv="results/summary/per_sample/{sample}.tsv",
    params:
        delimiter="\t",
        sample_name=lambda wildcards: wildcards.sample,
    conda:
        "../envs/python.yaml"
    log:
        "logs/summary/summary_results/{sample}.log",
    script:
        "../scripts/summary.py"
