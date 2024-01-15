
rule qualimap__mapping_quality_report:
    input:
        bam="results/self_contamination/{sample}/markdup.bam",
        bai="results/self_contamination/{sample}/markdup.bam.bai",
    output:
        report_dir=report(
            directory("results/self_contamination/{sample}/markdup/bamqc"),
            category="{sample}",
            labels={
                "Type": "Qualimap for markdup",
            },
            htmlindex="qualimapReport.html",
        ),
    params:
        extra=[
            "--paint-chromosome-limits",
            "-outformat PDF:HTML",
        ],
    resources:
        mem_mb=get_mem_mb_for_mapping_postprocess,
    log:
        "logs/qualimap/mapping_quality_report/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.7/wrappers/qualimap/bamqc"


rule check_coverage_from_qualimap:
    input:
        "results/self_contamination/{sample}/markdup/bamqc",
    output:
        "results/checks/{sample}/coverage_check.txt",
    params:
        genome_results_file=lambda wildcards, input: os.path.join(input[0], "genome_results.txt"),
        warn_threshold=config["coverage_check"]["warn_threshold"],
        fail_threshold=config["coverage_check"]["fail_threshold"],
    log:
        "logs/checks/coverage_check/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/coverage_check.py"
