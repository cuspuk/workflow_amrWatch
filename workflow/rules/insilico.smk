
rule prepare_tsv_insilico:
    output:
        "results/.primers.tsv",
    conda:
        "../envs/coreutils.yaml"
    params:
        primers="\n".join(["\t".join(t.split(",")) for t in config["in_silico_PCR"]["primers"]]),
    log:
        "logs/custom/prepare_tsv_insilico.log",
    localrule: True
    shell:
        "echo {params.primers:q} > {output} 2> {log}"


rule in_silico_PCR:
    input:
        fasta=infer_assembly_fasta,
        primers="results/.primers.tsv",
    output:
        "results/in_silico_PCR/{sample}.bed",
    params:
        bed="--bed" if config["in_silico_PCR"]["output_amplicon"] else "",
        mismatches=(
            f'-m {config["in_silico_PCR"]["max_mismatches"]}'
            if config["in_silico_PCR"].get("max_mismatches", None) is not None
            else ""
        ),
    conda:
        "../envs/seqkit_older.yaml"
    log:
        "logs/custom/in_silico_PCR/{sample}.log",
    shell:
        "(cat {input.fasta} | seqkit amplicon {params.mismatches} -p {input.primers} {params.bed} > {output}) 2> {log}"
