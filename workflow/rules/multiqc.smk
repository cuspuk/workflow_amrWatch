
rule multiqc__report:
    input:
        cutadapt=expand("results/reads/trimmed/{sample}.qc.txt", sample=get_sample_names()),
        fastqc=expand(
            "results/reads/trimmed/fastqc/{sample}_{pair}/fastqc_data.txt",
            sample=get_sample_names(),
            pair=["R1", "R2"],
        ),
        kraken=expand("results/kraken/{sample}.kreport2", sample=get_sample_names()),
    output:
        "results/summary/multiqc.html",
    params:
        use_input_files_only=True,
        extra=f"--config {workflow.basedir}/resources/multiqc.yaml",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v3.3.0/bio/multiqc"
