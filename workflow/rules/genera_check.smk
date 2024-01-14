rule curl__download_kraken_db:
    output:
        protected("{prefix_dir}/hash.k2d"),
    params:
        url=lambda wildcards, output: "https://genome-idx.s3.amazonaws.com/kraken/{tag}.tar.gz".format(
            tag=os.path.basename(os.path.dirname(output[0]))
        ),
        dirpath=lambda wildcards, output: os.path.dirname(output[0]),
    retries: 1
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
        "results/kraken/{sample}.bracken",
    params:
        db_dir=lambda wildcards, input: os.path.dirname(input.kraken_tax),
        read_length=config["foreign_contamination"]["read_length"],
        threshold=config["foreign_contamination"]["threshold"],
        classification_level=config["foreign_contamination"]["classification_level"],
    threads: min(config["threads"]["bracken"], config["max_threads"])
    log:
        "logs/kraken/analysis/{sample}.log",
    conda:
        "../envs/bracken.yaml"
    shell:
        "bracken -d {params.db_dir} -i {input.report} -o {output} -r {params.read_length}"
        " -t {params.threshold} -l {params.classification_level} > {log} 2>&1"


rule check_foreign_contamination:
    input:
        "results/kraken/{sample}.bracken",
    output:
        "results/checks/{sample}/foreign_contamination.txt",
    params:
        fraction_threshold=0.01,
    log:
        "logs/checks/foreign_contamination/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/genera_check.py"
