
rule multiqc__report:
    input:
        cutadapt="results/reads/trimmed/{sample}_cutadapt.json",
        fastqc=expand("results/reads/trimmed/fastqc/{{sample}}_{pair}.html", pair=["R1", "R2"]),
        bracken="results/kraken/{sample}.bracken",
    output:
        "results/multiqc/{sample}.html",
    params:
        extra="",  # Optional: extra parameters for multiqc.
        use_input_files_only=True,
    log:
        "logs/multiqc/{sample}.log",
    wrapper:
        "v3.3.0/bio/multiqc"
