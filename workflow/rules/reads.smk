rule cutadapt__trim:
    input:
        infer_fastqs_for_trimming,
    output:
        fastq1="results/reads/trimmed/{sample}_R1.fastq.gz",
        fastq2="results/reads/trimmed/{sample}_R2.fastq.gz",
        qc="results/reads/trimmed/{sample}.qc.txt",
    params:
        extra=get_cutadapt_extra_pe(),
    resources:
        mem_mb=get_mem_mb_for_trimming,
    threads: min(config["threads"]["reads__trimming"], config["max_threads"])
    log:
        "logs/cutadapt/trim_reads_pe/{sample}.log",
    wrapper:
        "v3.4.1/bio/cutadapt/pe"


rule fastqc__report:
    input:
        read=infer_fastq_path_for_fastqc,
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
        "https://github.com/xsitarcik/wrappers/raw/v1.12.12/wrappers/fastqc/quality"


rule kraken__download_db:
    output:
        protected("{prefix_dir}/hash.k2d"),
    params:
        url=lambda wildcards, output: "https://genome-idx.s3.amazonaws.com/kraken/{tag}.tar.gz".format(
            tag=os.path.basename(os.path.dirname(output[0]))
        ),
        dirpath=lambda wildcards, output: os.path.dirname(output[0]),
    retries: 1
    localrule: True
    log:
        "{prefix_dir}/logs/download.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "(mkdir -p {params.dirpath} && curl -SL {params.url} | tar zxvf - -C {params.dirpath}) > {log} 2>&1"


rule kraken__analysis:
    input:
        kraken_tax=os.path.join(config["foreign_contamination"]["kraken_dir"], "hash.k2d"),
        r1="results/reads/trimmed/{sample}_R1.fastq.gz",
        r2="results/reads/trimmed/{sample}_R2.fastq.gz",
    output:
        kraken_output=temp("results/kraken/{sample}.kraken"),
        report="results/kraken/{sample}.kreport2",
    params:
        save_memory="--memory-mapping" if config["foreign_contamination"]["save_memory"] else "",
        db_dir=lambda wildcards, input: os.path.dirname(input.kraken_tax),
    threads: min(config["threads"]["kraken"], config["max_threads"])
    log:
        "logs/kraken/analysis/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "(kraken2 --db {params.db_dir} --threads {threads} --paired --gzip-compressed"
        " {params.save_memory} --report {output.report} {input.r1} {input.r2} 1> {output.kraken_output}) 2> {log}"


rule bracken__correction:
    input:
        kraken_tax=os.path.join(config["foreign_contamination"]["kraken_dir"], "hash.k2d"),
        report="results/kraken/{sample}.kreport2",
    output:
        bracken="results/kraken/{sample}.bracken",
        kreport="results/kraken/{sample}_bracken_genuses.kreport2",
    params:
        db_dir=lambda wildcards, input: os.path.dirname(input.kraken_tax),
        read_length=config["foreign_contamination"]["read_length"],
        threshold=config["foreign_contamination"]["bracken_threshold"],
        classification_level=config["foreign_contamination"]["classification_level"],
    threads: min(config["threads"]["bracken"], config["max_threads"])
    log:
        "logs/kraken/analysis/{sample}.log",
    conda:
        "../envs/bracken.yaml"
    shell:
        "bracken -d {params.db_dir} -i {input.report} -o {output.bracken} -r {params.read_length}"
        " -t {params.threshold} -l {params.classification_level} > {log} 2>&1"


rule multiqc__report:
    input:
        cutadapt=expand("results/reads/trimmed/{sample}.qc.txt", sample=get_sample_names_with_reads_as_input()),
        fastqc=expand(
            "results/reads/trimmed/fastqc/{sample}_{pair}/fastqc_data.txt",
            sample=get_sample_names_with_reads_as_input(),
            pair=["R1", "R2"],
        ),
        bracken=expand(
            "results/kraken/{sample}_bracken_genuses.kreport2", sample=get_sample_names_with_reads_as_input()
        ),
    output:
        "results/summary/multiqc.html",
    params:
        use_input_files_only=True,
        extra=f"--config {workflow.basedir}/resources/multiqc.yaml",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v3.4.1/bio/multiqc"
