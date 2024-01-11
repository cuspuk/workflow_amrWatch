checkpoint assembly_constructed:
    input:
        "results/assembly/{sample}/assembly.gfa",
    output:
        "results/checks/{sample}/assembly.txt",
    conda:
        "../envs/coreutils.yaml"
    log:
        "logs/checks/assembly/{sample}.log",
    shell:
        "(([ -s {input} ] && echo 'PASS') || echo 'FAIL') > {output} 2> {log}"


rule request_bandage_if_assembly_exists:
    input:
        optional_bandage_outputs,
    output:
        "results/checks/{sample}/.bandage_requested.txt",
    conda:
        "../envs/coreutils.yaml"
    log:
        "logs/checks/assembly/{sample}.log",
    shell:
        "touch {output} 2> {log}"


# rule aggregate:
#     input:
#         expand("results/checks/{sample}/assembly.txt",sample=)
#     output:
#         "results/.aggregation/{sample}.txt"
#     shell:
#         "touch {output}"


rule summary_all_checks:
    input:
        foreign_contamination="results/checks/{sample}/foreign_contamination.txt",
        assembly_check="results/checks/{sample}/assembly.txt",
    output:
        "results/checks/{sample}/summary.txt",
    log:
        "logs/checks/summary/{sample}.log",
    conda:
        "../envs/bracken.yaml"
    shell:
        "echo PASS > {output}"  # > {log} 2>&1"
