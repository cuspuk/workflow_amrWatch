rule bwa__build_index:
    input:
        "results/assembly/{sample}/assembly_cleaned.fasta",
    output:
        idx=temp(multiext("results/assembly/{sample}/bwa_index/bwa", ".amb", ".ann", ".bwt", ".pac", ".sa")),
    params:
        extra="-a bwtsw",
    log:
        "logs/self_contamination/bwa_index/{sample}.log",
    wrapper:
        "v3.4.1/bio/bwa/index"


rule bwa__map_to_assembly:
    input:
        reads=["results/reads/trimmed/{sample}_R1.fastq.gz", "results/reads/trimmed/{sample}_R2.fastq.gz"],
        idx=multiext("results/assembly/{sample}/bwa_index/bwa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("results/self_contamination/{sample}/mapped.bam"),
    log:
        "logs/self_contamination/bwa/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",
        sort_order="queryname",
        sort_extra="",
    threads: min(config["threads"]["mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    wrapper:
        "v3.4.1/bio/bwa/mem"


rule samtools__filter_unmapped:
    input:
        "results/self_contamination/{sample}/mapped.bam",
    output:
        bam=temp("results/self_contamination/{sample}/filtered.bam"),
    log:
        "logs/self_contamination/samtools_filter/{sample}.log",
    params:
        extra="-F '0x4'",
        region="",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    wrapper:
        "v3.4.1/bio/samtools/view"


rule samtools__fixmate:
    input:
        "results/self_contamination/{sample}/mapped.bam",
    output:
        temp("results/self_contamination/{sample}/fixmate.bam"),
    params:
        extra="-m",
    log:
        "logs/self_contamination/samtools_fixmate/{sample}.log",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    wrapper:
        "v3.4.1/bio/samtools/fixmate/"


rule samtools__sort_after_fixmate:
    input:
        "results/self_contamination/{sample}/fixmate.bam",
    output:
        temp("results/self_contamination/{sample}/sorted.bam"),
    log:
        "logs/self_contamination/samtools_sort/{sample}.log",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    wrapper:
        "v3.4.1/bio/samtools/sort"


rule samtools__markdup:
    input:
        "results/self_contamination/{sample}/sorted.bam",
    output:
        temp("results/self_contamination/{sample}/markdup.bam"),
    log:
        "logs/self_contamination/samtools_markdup/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    shell:
        "samtools markdup -@ {threads} {input} {output} > {log} 2>&1"


rule samtools__index:
    input:
        "results/self_contamination/{sample}/markdup.bam",
    output:
        temp("results/self_contamination/{sample}/markdup.bam.bai"),
    log:
        "logs/self_contamination/samtools_index/{sample}.log",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    wrapper:
        "v3.4.1/bio/samtools/index"


rule qualimap__report:
    input:
        bam="results/self_contamination/{sample}/markdup.bam",
        bai="results/self_contamination/{sample}/markdup.bam.bai",
    output:
        report_dir=directory("results/self_contamination/{sample}/markdup/bamqc"),
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    log:
        "logs/qualimap/mapping_quality_report/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.12/wrappers/qualimap/bamqc"


rule samtools__faidx:
    input:
        "results/assembly/{sample}/assembly_cleaned.fasta",
    output:
        "results/assembly/{sample}/assembly_cleaned.fasta.fai",
    log:
        "logs/self_contamination/samtools_faidx/{sample}.log",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    wrapper:
        "v3.4.1/bio/samtools/faidx"


rule bcftools__mpileup:
    input:
        alignments="results/self_contamination/{sample}/markdup.bam",
        ref="results/assembly/{sample}/assembly_cleaned.fasta",
        index="results/assembly/{sample}/assembly_cleaned.fasta.fai",
    output:
        pileup=temp("results/self_contamination/{sample}/piled.bcf"),
    params:
        uncompressed_bcf=True,
        extra="--gvcf 20 --no-BAQ -d 1000 --min-MQ 60",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    log:
        "logs/self_contamination/bcftools_mpileup/{sample}.log",
    wrapper:
        "v3.4.1/bio/bcftools/mpileup"


rule bcftools__variants:
    input:
        pileup="results/self_contamination/{sample}/piled.bcf",
    output:
        calls=temp("results/self_contamination/{sample}/called.bcf"),
    params:
        uncompressed_bcf=True,
        caller="--multiallelic-caller",
        extra="-v",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    log:
        "logs/self_contamination/bcftools_call/{sample}.log",
    wrapper:
        "v3.4.1/bio/bcftools/call"


rule bcftools__norm:
    input:
        "results/self_contamination/{sample}/called.bcf",
    output:
        temp("results/self_contamination/{sample}/norm.bcf"),
    params:
        uncompressed_bcf=True,
        extra="--rm-dup all",
    log:
        "logs/self_contamination/bcftools_norm/{sample}.log",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    wrapper:
        "v3.4.1/bio/bcftools/norm"


rule bcftools__filter:
    input:
        "results/self_contamination/{sample}/norm.bcf",
    output:
        "results/self_contamination/{sample}/filtered.vcf",
    log:
        "logs/self_contamination/bcftools_filter/{sample}.log",
    params:
        filter="-i 'MAF > {polymorph_rate}'".format(polymorph_rate=config["self_contamination"]["polymorph_rate"]),
        extra="",
    threads: min(config["threads"]["mapping_postprocess"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    wrapper:
        "v3.4.1/bio/bcftools/filter"
