
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
