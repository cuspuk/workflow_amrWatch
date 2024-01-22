rule gtdbtk__classify:
    input:
        assembly=infer_assembly_fasta,
        gtdb=os.path.join(config["gtdb_dirpath"], "db"),
    output:
        dir=directory("results/taxonomy/{sample}/classify"),
    params:
        assembly_dir=lambda wildcards, input: os.path.dirname(input.assembly),
        out_dir=lambda wildcards, output: os.path.dirname(output.dir),
    threads: min(config["threads"]["gtdb__classify"], config["max_threads"])
    conda:
        "../envs/gtdbtk.yaml"
    resources:
        mem_mb=get_mem_mb_for_gtdb,
    log:
        "logs/taxonomy/gtdb_classify/{sample}.log",
    shell:
        "(export GTDBTK_DATA_PATH={input.gtdb:q} "
        " && echo '{input.assembly}\t{wildcards.sample}' > $TMPDIR/batchfile_{wildcards.sample}.txt"
        " && gtdbtk classify_wf --batchfile $TMPDIR/batchfile_{wildcards.sample}.txt"
        " --cpus {threads} --tmpdir $TMPDIR --extension fasta --out_dir {params.out_dir} --mash_db {input.gtdb}) > {log} 2>&1"


checkpoint gtdbtk__parse_taxa:
    input:
        gtdb_outdir="results/taxonomy/{sample}/classify",
    output:
        "results/taxonomy/{sample}/parsed_taxa.txt",
    params:
        gtdb_tsv=lambda wildcards, input: os.path.join(input.gtdb_outdir, "gtdbtk.bac120.summary.tsv"),
    conda:
        "../envs/coreutils.yaml"
    localrule: True
    log:
        "logs/taxonomy/parse_taxa/{sample}.log",
    shell:
        '(cut -f2 {params.gtdb_tsv} | tail -n 1 | sed -e "s/.*;s__//") > {output} 2> {log}'
