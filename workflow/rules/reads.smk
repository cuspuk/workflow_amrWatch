rule cutadapt__trim_reads_pe:
    input:
        get_reads_for_trimming,
    output:
        fastq1=temp("results/reads/trimmed/{sample}_R1.fastq.gz"),
        fastq2=temp("results/reads/trimmed/{sample}_R2.fastq.gz"),
        qc=temp("results/reads/trimmed/{sample}.qc.txt"),
    params:
        extra=get_cutadapt_extra_pe(),
    resources:
        mem_mb=get_mem_mb_for_trimming,
    threads: min(config["threads"]["reads__trimming"], config["max_threads"])
    log:
        "logs/cutadapt/trim_reads_pe/{sample}.log",
    wrapper:
        "v3.3.3/bio/cutadapt/pe"


rule fastqc__quality_report:
    input:
        read=infer_fastq_path,
    output:
        html=report(
            "results/reads/{step}/fastqc/{sample}_{pair}.html",
            category="{sample}",
            labels={
                "Type": "Fastqc {pair} - {step}",
            },
        ),
        zip="results/reads/{step}/fastqc/{sample}_{pair}.zip",
        qc_data="results/reads/{step}/fastqc/{sample}_{pair}/fastqc_data.txt",
        summary_txt="results/reads/{step}/fastqc/{sample}_{pair}/summary.txt",
    threads: min(config["threads"]["fastqc"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_fastqc,
    log:
        "logs/fastqc/{step}/{sample}_{pair}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.6/wrappers/fastqc/quality"
