rule amrfinder__download_db:
    output:
        db=directory(config["amrfinder"]["db_dir"]),
    params:
        db_parent=lambda wildcards, output: os.path.dirname(output.db),
    conda:
        "../envs/amrfinder.yaml"
    log:
        os.path.join(os.path.dirname(config["amrfinder"]["db_dir"]), "logs", "download.log"),
    shell:
        "(mkdir -p {params.db_parent} && amrfinder_update -d {params.db_parent}) > {log} 2>&1"


rule amrfinder__call:
    input:
        contigs=infer_assembly_fasta,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
        db=config["amrfinder"]["db_dir"],
    output:
        tsv="results/amr_detect/{sample}/amrfinder.tsv",
    params:
        organism_arg=get_organism_for_amrfinder,
        plus_arg="--plus" if config["amrfinder"]["plus"] else "",
    threads: min(config["threads"]["amrfinder"], config["max_threads"])
    conda:
        "../envs/amrfinder.yaml"
    log:
        "logs/amr_detect/amrfinder/{sample}.log",
    shell:
        "amrfinder {params.plus_arg} -d {input.db} --nucleotide {input.contigs}"
        " --threads {threads} {params.organism_arg} -o {output.tsv} > {log} 2>&1"


rule mlst__call:
    input:
        contigs=infer_assembly_fasta,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
        blast_mlst=os.path.join(config["mlst_db_dir"], "blast", "mlst.fa"),
        pubmlst=os.path.join(config["mlst_db_dir"], "pubmlst"),
    output:
        "results/amr_detect/{sample}/mlst.tsv",
    params:
        scheme_arg=get_taxonomy_for_mlst,
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/amr_detect/mlst/{sample}.log",
    shell:
        "mlst {input.contigs} --blastdb {input.blast_mlst} --datadir {input.pubmlst} {params.scheme_arg} --legacy > {output} 2> {log}"


rule abricate_download_db:
    output:
        db=protected(os.path.join(config["abricate"]["db_dir"], "sequences")),
    params:
        db_name=lambda wildcards, output: os.path.basename(os.path.dirname(output.db)),
        db_dir=lambda wildcards, output: os.path.dirname(os.path.dirname(output.db)),
    localrule: True
    conda:
        "../envs/abricate.yaml"
    log:
        os.path.join(config["abricate"]["db_dir"], "logs", "download.log"),
    shell:
        "abricate-get_db --db={params.db_name} --dbdir={params.db_dir} > {log} 2>&1"


rule abricate__call:
    input:
        fasta=infer_assembly_fasta,
        db=os.path.join(config["abricate"]["db_dir"], "sequences"),
    output:
        "results/amr_detect/{sample}/abricate.tsv",
    params:
        db_name=lambda wildcards, input: os.path.basename(os.path.dirname(input.db)),
        db_dir=lambda wildcards, input: os.path.dirname(os.path.dirname(input.db)),
        min_identity=config["abricate"]["min_identity"],
        min_coverage=config["abricate"]["min_coverage"],
    conda:
        "../envs/abricate.yaml"
    log:
        "logs/amr_detect/abricate/{sample}.log",
    shell:
        "abricate --datadir {params.db_dir} --db {params.db_name} --minid {params.min_identity} --mincov {params.min_coverage}"
        " {input.fasta} > {output} 2> {log}"


rule kleborate__call:
    input:
        fasta=infer_assembly_fasta,
    output:
        tsv="results/amr_detect/{sample}/kleborate/klebsiella_pneumo_complex_output.txt",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.tsv),
    conda:
        "../envs/kleborate.yaml"
    log:
        "logs/amr_detect/kleborate/{sample}.log",
    shell:
        "kleborate -o {params.outdir} -a {input.fasta} -p kpsc --trim_headers > {log} 2>&1"


rule spatyper__download_db:
    output:
        repeats=protected(os.path.join(config["spatyper_db_dir"], "sparepeats.fasta")),
        order=protected(os.path.join(config["spatyper_db_dir"], "spatypes.txt")),
    params:
        db_dir=lambda wildcards, output: os.path.dirname(output.repeats),
    localrule: True
    conda:
        "../envs/spatyper.yaml"
    log:
        os.path.join(config["spatyper_db_dir"], "logs", "download.log"),
    script:
        "../scripts/spatyper_db_download.py"


rule spatyper__call:
    input:
        fasta=infer_assembly_fasta,
        repeats=os.path.join(config["spatyper_db_dir"], "sparepeats.fasta"),
        order=os.path.join(config["spatyper_db_dir"], "spatypes.txt"),
    output:
        "results/amr_detect/{sample}/spa_typer.tsv",
    conda:
        "../envs/spatyper.yaml"
    log:
        "logs/amr_detect/spa_typer/{sample}.log",
    shell:
        "spaTyper -f {input.fasta} -r {input.repeats} -o {input.order} --output {output} > {log} 2>&1"


rule etoki__call:
    input:
        infer_assembly_fasta,
    output:
        "results/amr_detect/{sample}/etoki_ebeis.tsv",
    conda:
        "../envs/etoki.yaml"
    log:
        "logs/amr_detect/etoki_ebeis/{sample}.log",
    shell:
        "EToKi.py EBEis -q {input} > {output} 2>{log}"


rule etoki__iscrispol:
    input:
        infer_assembly_fasta,
    output:
        "results/amr_detect/{sample}/crispol.tsv",
    conda:
        "../envs/etoki.yaml"
    log:
        "logs/amr_detect/etoki_iscrispol/{sample}.log",
    shell:
        "EToKi.py isCRISPOL {input} > {output} 2>{log}"


rule sccmec__download_db:
    output:
        db_proteins=protected(os.path.join(config["SCCmec_db_dir"], "proteins.fasta")),
        db_subtypes=protected(os.path.join(config["SCCmec_db_dir"], "subtypes.fasta")),
        db_primers=protected(os.path.join(config["SCCmec_db_dir"], "primers.fasta")),
    params:
        repo="https://github.com/staphopia/staphopia-sccmec/archive/refs/heads/master.zip",
        repo_db_name="staphopia-sccmec-master/share/staphopia-sccmec/data",
        db_dir=lambda wildcards, output: os.path.dirname(output.db_proteins),
    conda:
        "../envs/curl_with_unzip.yaml"
    log:
        os.path.join(config["SCCmec_db_dir"], "logs", "download.log"),
    script:
        "../scripts/sccmec_download_db.py"


rule sccmec__call:
    input:
        assembly=infer_assembly_fasta,
        db_proteins=os.path.join(config["SCCmec_db_dir"], "proteins.fasta"),
        db_subtypes=os.path.join(config["SCCmec_db_dir"], "subtypes.fasta"),
        db_primers=os.path.join(config["SCCmec_db_dir"], "primers.fasta"),
    output:
        "results/amr_detect/{sample}/SCCmec.tsv",
    params:
        ext=lambda wildcards, input: os.path.splitext(input[0])[1],
        db_dir=lambda wildcards, input: os.path.dirname(input.db_proteins),
    conda:
        "../envs/sccmec.yaml"
    log:
        "logs/amr_detect/{sample}/sccmec.log",
    shell:
        "staphopia-sccmec --assembly {input.assembly} --sccmec {params.db_dir} --ext {params.ext} > {output} 2>{log}"


rule mob_suite__download_db:
    output:
        db=protected(os.path.join(config["ncbi_plasmid_db_dir"], "ncbi_plasmid_full_seqs.fas.msh")),
    params:
        db_dir=lambda wildcards, output: os.path.dirname(output.db),
    conda:
        "../envs/mob_suite.yaml"
    log:
        "logs/plasmids/download_db.log",
    shell:
        "mob_init --database_directory {params.db_dir} > {log} 2>&1"


rule mob_suite__typer:
    input:
        fasta=infer_assembly_fasta,
        db=os.path.join(config["ncbi_plasmid_db_dir"], "ncbi_plasmid_full_seqs.fas.msh"),
    output:
        "results/plasmids/{sample}/mob_typer.txt",
    params:
        db_dir=lambda wildcards, input: os.path.dirname(input.db),
    conda:
        "../envs/mob_suite.yaml"
    log:
        "logs/plasmids/{sample}/mob_typer.log",
    shell:
        "mob_typer --infile {input.fasta} --database_directory {params.db_dir} --out_file {output} > {log} 2>&1"


rule sistr__cmd__call:
    input:
        infer_assembly_fasta,
    output:
        serovar="results/amr_detect/{sample}/sistr_serovar.tab",
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.serovar),
    conda:
        "../envs/sistr_cmd.yaml"
    log:
        "logs/amr_detect/sistr/{sample}.log",
    shell:
        "(mkdir -p {params.out_dir} && sistr --qc --output-format tab --output-prediction {output.serovar}"
        " --no-cgmlst {input}) > {log} 2>&1"


rule rgi__download_db:
    output:
        json=protected(os.path.join(config["rgi_db_dir"], "card.json")),
    params:
        db_dir=lambda wildcards, output: os.path.dirname(output.json),
        db_url="https://card.mcmaster.ca/latest/data",
    conda:
        "../envs/curl.yaml"
    localrule: True
    log:
        os.path.join(config["rgi_db_dir"], "logs", "download.log"),
    shell:
        "mkdir -p {params.db_dir} && (curl -SL {params.db_url} | tar xjvf - -C {params.db_dir}) > {log} 2>&1"


rule rgi__load_db:
    input:
        json=os.path.join(config["rgi_db_dir"], "card.json"),
    output:
        loaded_db=directory("localDB"),
    conda:
        "../envs/rgi.yaml"
    log:
        "logs/amr_detect/rgi_db_load.log",
    shell:
        "(rm -rf {output.loaded_db} && rgi load --card_json {input.json} --local) > {log} 2>&1"


rule rgi__call:
    input:
        json=os.path.join(config["rgi_db_dir"], "card.json"),
        loaded_db="localDB",
        assembly=infer_assembly_fasta,
    output:
        txt="results/amr_detect/{sample}/rgi_main.txt",
        json="results/amr_detect/{sample}/rgi_main.json",
    params:
        extra="--include_nudge",
        out_prefix=lambda wildcards, output: os.path.splitext(output.txt)[0],
    conda:
        "../envs/rgi.yaml"
    threads: min(config["threads"]["amrfinder"], config["max_threads"])
    log:
        "logs/amr_detect/rgi/{sample}.log",
    shell:
        "rgi main --input_sequence {input.assembly} --local --clean {params.extra}"
        " --output_file {params.out_prefix} --num_threads {threads} --split_prodigal_jobs > {log} 2>&1"


rule resfinder__download_db:
    output:
        resfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "resfinder_db", "VERSION"),
        pointfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "pointfinder_db", "VERSION"),
        disinfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "disinfinder_db", "VERSION"),
    params:
        resfinder_db_url="https://bitbucket.org/genomicepidemiology/resfinder_db/",
        pointfinder_db_url="https://bitbucket.org/genomicepidemiology/pointfinder_db/",
        disinfinder_db_url="https://bitbucket.org/genomicepidemiology/disinfinder_db/",
        resfinder_db_dir=lambda wildcards, output: os.path.dirname(output.resfinder_db_ver),
        pointfinder_db_dir=lambda wildcards, output: os.path.dirname(output.pointfinder_db_ver),
        disinfinder_db_dir=lambda wildcards, output: os.path.dirname(output.disinfinder_db_ver),
    conda:
        "../envs/git.yaml"
    localrule: True
    log:
        os.path.join(config["resfinder"]["db_dir"], "logs", "download.log"),
    shell:
        "( git clone {params.resfinder_db_url} {params.resfinder_db_dir}"
        " && git clone {params.pointfinder_db_url} {params.pointfinder_db_dir}"
        " && git clone {params.disinfinder_db_url} {params.disinfinder_db_dir}"
        " ) > {log} 2>&1"


rule resfinder__kma_index:
    input:
        resfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "resfinder_db", "VERSION"),
        pointfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "pointfinder_db", "VERSION"),
        disinfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "disinfinder_db", "VERSION"),
    output:
        resfinder_out=[
            os.path.join(config["resfinder"]["db_dir"], "resfinder_db", "{value}.comp.b").format(value=value)
            for value in [
                "misc",
                "pseudomonicacid",
                "fusidicacid",
                "phenicol",
                "glycopeptide",
                "trimethoprim",
                "oxazolidinone",
                "tetracycline",
                "quinolone",
                "nitroimidazole",
                "fosfomycin",
                "aminoglycoside",
                "macrolide",
                "sulphonamide",
                "rifampicin",
                "colistin",
                "beta-lactam",
            ]
        ],
        pointfinder_out=[
            os.path.join(config["resfinder"]["db_dir"], "pointfinder_db", "{value}", "{value}.comp.b").format(
                value=value
            )
            for value in [
                "campylobacter",
                "escherichia_coli",
                "enterococcus_faecalis",
                "enterococcus_faecium",
                "neisseria_gonorrhoeae",
                "salmonella",
                "helicobacter_pylori",
                "klebsiella",
                "plasmodium_falciparum",
                "staphylococcus_aureus",
                "mycobacterium_tuberculosis",
            ]
        ],
        disinfinder_out=[os.path.join(config["resfinder"]["db_dir"], "disinfinder_db", "disinfectants.comp.b")],
    params:
        suffix=".comp.b",
        resfinder_db_dir=lambda wildcards, input: os.path.dirname(input.resfinder_db_ver),
        pointfinder_db_dir=lambda wildcards, input: os.path.dirname(input.pointfinder_db_ver),
        disinfinder_db_dir=lambda wildcards, input: os.path.dirname(input.disinfinder_db_ver),
    log:
        os.path.join(config["resfinder"]["db_dir"], "logs", "kma_index.log"),
    conda:
        "../envs/resfinder.yaml"
    script:
        "../scripts/kma_index.py"


rule resfinder__call:
    input:
        resfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "resfinder_db", "VERSION"),
        pointfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "pointfinder_db", "VERSION"),
        disinfinder_db_ver=os.path.join(config["resfinder"]["db_dir"], "disinfinder_db", "VERSION"),
        inferred_input=infer_resfinder_input,
        taxa="results/taxonomy/{sample}/parsed_taxa.txt",
        kma_resfinder=expand(
            os.path.join(config["resfinder"]["db_dir"], "resfinder_db", "{value}.comp.b"),
            value=[
                "misc",
                "pseudomonicacid",
                "fusidicacid",
                "phenicol",
                "glycopeptide",
                "trimethoprim",
                "oxazolidinone",
                "tetracycline",
                "quinolone",
                "nitroimidazole",
                "fosfomycin",
                "aminoglycoside",
                "macrolide",
                "sulphonamide",
                "rifampicin",
                "colistin",
                "beta-lactam",
            ],
        ),
    output:
        tsv="results/amr_detect/{sample}/resfinder/ResFinder_results_tab.txt",
        pointfinder="results/amr_detect/{sample}/resfinder/PointFinder_results.txt",
    params:
        resfinder_db_dir=lambda wildcards, input: os.path.dirname(input.resfinder_db_ver),
        pointfinder_db_dir=lambda wildcards, input: os.path.dirname(input.pointfinder_db_ver),
        disinfinder_db_dir=lambda wildcards, input: os.path.dirname(input.disinfinder_db_ver),
        species=get_taxonomy_for_resfinder,
        outdir=lambda wildcards, output: os.path.dirname(output.tsv),
        input_arg=lambda wildcards, input: "--inputfasta" if isinstance(input.inferred_input, str) else "--inputfastq",
        min_cov=config["resfinder"]["min_coverage"],
        threshold=config["resfinder"]["threshold"],
    conda:
        "../envs/resfinder.yaml"
    log:
        "logs/amr_detect/resfinder/{sample}.log",
    shell:
        "(python -m resfinder {params.input_arg} {input.inferred_input} --ignore_missing_species"
        " --min_cov {params.min_cov} --threshold {params.threshold} -s {params.species:q}"
        " --disinfectant --db_path_disinf {params.disinfinder_db_dir}"
        " --point --db_path_point {params.pointfinder_db_dir}"
        " --acquired --db_path_res {params.resfinder_db_dir}"
        " -o {params.outdir} && touch {output.pointfinder}) > {log} 2>&1"


rule seqsero2__call:
    input:
        infer_assembly_fasta,
    output:
        tsv="results/amr_detect/{sample}/seqsero_summary.tsv",
    params:
        out_dir=lambda wildcards, output: os.path.join(os.path.dirname(output.tsv), "seqsero"),
        header="\t".join(
            [
                "Sample name",
                "Output directory",
                "Input files",
                "O antigen prediction",
                "H1 antigen prediction(fliC)",
                "H2 antigen prediction(fljB)",
                "Predicted identification",
                "Predicted antigenic profile",
                "Predicted serotype",
                "Note",
            ]
        ),
    conda:
        "../envs/seqsero.yaml"
    log:
        "logs/amr_detect/seqsero/{sample}.log",
    shell:
        "(mkdir -p {params.out_dir}"
        " && SeqSero2_package.py -t 4 -m k -s -i {input} -n {wildcards.sample} -d {params.out_dir}"
        " && echo -e {params.header:q} > {output.tsv}"
        " && cat {params.out_dir}/SeqSero_result.tsv >> {output.tsv}"
        " ) > {log} 2>&1"


rule clonal_complex_profiles_download:
    output:
        tsv=os.path.join(config["clonal_complex"]["db_dir"], "{ncbi_name}.tsv"),
    params:
        url=lambda wildcards: f"https://rest.pubmlst.org/db/pubmlst_{wildcards.ncbi_name}_seqdef/schemes/1/profiles_csv",
        out_dir=lambda wildcards, output: os.path.dirname(output.tsv),
    log:
        os.path.join(config["clonal_complex"]["db_dir"], "logs", "{ncbi_name}.log"),
    conda:
        "../envs/curl.yaml"
    shell:
        "(mkdir -p {params.out_dir} && curl -s -o {output.tsv} {params.url}) > {log} 2>&1"


rule mlst_clonal_complex:
    input:
        mlst="results/amr_detect/{sample}/mlst.tsv",
        profile_tsv=infer_profile_for_clonal_complex,
    output:
        "results/amr_detect/{sample}/clonal_complex.tsv",
    params:
        mlst_index=2,
    conda:
        "../envs/python.yaml"
    localrule: True
    log:
        "logs/amr_detect/mlst_clonal_complex/{sample}.log",
    script:
        "../scripts/mlst_to_clonal_complex.py"


rule pneumokity_download:
    output:
        src=os.path.join(config["pneumokity_source_dir"], "pneumokity.py"),
    params:
        repo="https://github.com/CarmenSheppard/PneumoKITy/archive/refs/tags/v1.0.zip",
        name_after_unzip="PneumoKITy-1.0",
        dir=lambda wildcards, output: os.path.dirname(output.src),
    conda:
        "../envs/curl_with_unzip.yaml"
    localrule: True
    log:
        os.path.join(config["pneumokity_source_dir"], "logs", "download.log"),
    script:
        "../scripts/pneumokity_download.py"


rule pneumokity_call:
    input:
        fasta=infer_assembly_fasta,
        src=os.path.join(config["pneumokity_source_dir"], "pneumokity.py"),
    output:
        csv="results/amr_detect/{sample}/pneumo_capsular_typing/{sample}_result_data.csv",
    params:
        outdir=lambda wildcards, output: os.path.dirname(os.path.dirname(output.csv)),
    conda:
        "../envs/pneumokity.yaml"
    log:
        "logs/amr_detect/pneumokity/{sample}.log",
    shell:
        "python3 {input.src} pure -a {input.fasta} -t {threads} -o {params.outdir} -s {wildcards.sample} > {log} 2>&1"
