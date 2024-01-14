
rule multiqc__report:
    input:
        cutadapt="results/reads/trimmed/{sample}_cutadapt.json",
        fastqc=expand("results/reads/trimmed/fastqc/{{sample}}_{pair}/fastqc_data.txt", pair=["R1", "R2"]),
        kraken="results/kraken/{sample}.kreport2",
    output:
        "results/multiqc/{sample}.html",
    params:
        use_input_files_only=True,
        extra=f"--config {workflow.basedir}/resources/multiqc.yaml",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/{sample}.log",
    wrapper:
        "v3.3.0/bio/multiqc"
