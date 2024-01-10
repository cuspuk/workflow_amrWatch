rule summary_all_checks:
    input:
        foreign_contamination="results/checks/{sample}/foreign_contamination.txt",
    output:
        "results/checks/{sample}/summary.txt",
    log:
        "logs/checks/summary/{sample}.log",
    conda:
        "../envs/bracken.yaml"
    shell:
        "echo PASS > {output}"  # > {log} 2>&1"
