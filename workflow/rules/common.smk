from snakemake.utils import validate
import json
import re


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")


pepfile: config["pepfile"]


validate(pep.sample_table, "../schemas/samples.schema.yaml")


def get_sample_names() -> list[str]:
    return list(pep.sample_table["sample_name"].values)


def sample_has_asssembly_as_input(sample: str) -> bool:
    try:
        assembly = pep.sample_table.loc[sample][["fasta"]].iloc[0]
        return True if assembly else False
    except KeyError:
        return False


def sample_has_long_reads(sample: str) -> bool:
    try:
        long_reads = pep.sample_table.loc[sample][["long"]].iloc[0]
        return True if long_reads else False
    except KeyError:
        return False


def validate_dynamic_config():
    if config["resfinder"]["input_to_use"] == "reads":
        samples_with_assembly_as_input = [
            sample for sample in get_sample_names() if sample_has_asssembly_as_input(sample)
        ]
        if samples_with_assembly_as_input:
            raise ValueError(
                f"resfinder input_to_use is set to reads, but pepfile gives that input is assembly not reads. Relevant samples: {samples_with_assembly_as_input}"
            )


validate_dynamic_config()


def get_sample_names_with_reads_as_input() -> list[str]:
    return [sample for sample in get_sample_names() if not sample_has_asssembly_as_input(sample)]


def get_constraints() -> dict[str, list[str]]:
    return {"sample": "|".join(get_sample_names()), "pair": ["R1", "R2"]}


def get_first_fastq_for_sample_from_pep(sample: str, read_pair="fq1") -> str:
    return pep.sample_table.loc[sample][[read_pair]][0]


def get_fasta_for_sample_from_pep(sample: str) -> str:
    return pep.sample_table.loc[sample][["fasta"]][0]


def get_long_reads_for_sample_from_pep(sample: str) -> str:
    return pep.sample_table.loc[sample][["long"]][0]


with open(f"{workflow.basedir}/resources/gtdb_amrfinder.json", "r") as f:
    AMRFINDER_MAP = json.load(f)

with open(f"{workflow.basedir}/resources/gtdb_mlst.json", "r") as f:
    MLST_MAP = json.load(f)


def get_key_for_value_from_db(value: str, db: dict) -> str:
    # search species first
    for key in db:
        if "species" not in db[key]:
            continue
        pattern = f'{db[key]["genus"]} {db[key]["species"]}'
        if re.match(pattern, value):
            return key

    # search genus
    for key in db:
        if "species" in db[key]:
            continue
        if re.match(db[key]["genus"], value):
            return key
    raise KeyError


def get_parsed_taxa_from_gtdbtk_for_sample(sample: str):
    with checkpoints.checkpoint_parse_taxa_gtdbtk.get(sample=sample).output[0].open() as f:
        return f.read().strip()


def check_preassembly_QC_for_sample(sample: str) -> bool:
    with checkpoints.checkpoint_pre_assembly_QC.get(sample=sample).output[0].open() as f:
        return all([line.startswith(("PASS", "WARN")) for line in f.readlines()])


def check_assembly_construction_success_for_sample(sample: str):
    with checkpoints.checkpoint_assembly_construction.get(sample=sample).output[0].open() as f:
        return all([line.startswith(("PASS", "WARN")) for line in f.readlines()])


def check_all_checks_success_for_sample(sample: str):
    with checkpoints.checkpoint_request_post_assembly_checks_if_relevant.get(sample=sample).output[0].open() as f:
        header = f.readline()
        return all([line.startswith(("PASS", "WARN")) for line in f.readlines()])


def get_sample_names_passing_all_checks():
    sample_names = get_sample_names()
    return [s for s in sample_names if check_all_checks_success_for_sample(s)]


def infer_amr_detection_results_for_harmonize(wildcards):
    sample_names = get_sample_names_passing_all_checks()
    rgi = [f"results/hamronization/rgi/{sample}.tsv" for sample in sample_names]
    amrfinder = [f"results/hamronization/amrfinder/{sample}.tsv" for sample in sample_names]
    # abricate = [f"results/hamronization/abricate/{sample}.tsv" for sample in sample_names]
    if config["resfinder"]["input_to_use"] == "assembly":
        resfinder = [f"results/hamronization/resfinder/{sample}.tsv" for sample in sample_names]
        pointfinder = [f"results/hamronization/pointfinder/{sample}.tsv" for sample in sample_names]
    else:
        resfinder = []
        pointfinder = []
    return rgi + amrfinder + resfinder + pointfinder


def request_hamronize_or_nothing(wildcards):
    if not config["run_hamronization"]:
        return "results/hamronization/hamronization_skipped.txt"
    passed_samples = get_sample_names_passing_all_checks()
    if len(passed_samples) > 1:
        return expand("results/hamronization/summary.{ext}", ext=["tsv", "html"])
    else:
        return "results/hamronization/hamronization_skipped.txt"


def get_outputs():
    sample_names = get_sample_names()
    outputs = {
        "final_results": expand("results/checks/{sample}/.final_results_requested.tsv", sample=sample_names),
        "summary": "results/summary/summary.tsv",
        "amrfinder_summary": "results/summary/amrfinder_output.tsv",
    }

    if config["in_silico_PCR"]["do"]:
        outputs["in_silico_PCR"] = expand("results/in_silico_PCR/{sample}.bed", sample=sample_names)

    if len(sample_names) > 1:
        outputs["hamronization"] = "results/hamronization/hamronization_requested.txt"

    if samples_with_reads := get_sample_names_with_reads_as_input():
        if len(samples_with_reads) > 1:
            outputs["multiqc"] = "results/summary/multiqc.html"
        else:
            sample = samples_with_reads[0]
            outputs["qc"] = [
                f"results/reads/trimmed/fastqc/{sample}_R1/fastqc_data.txt",
                f"results/reads/trimmed/fastqc/{sample}_R2/fastqc_data.txt",
                f"results/kraken/{sample}.bracken",
            ]
    return outputs


def get_taxonomy_dependant_outputs(sample: str, taxa: str) -> dict[str, str]:
    outputs: dict[str, str] = {}
    if taxa.startswith("Klebsiella"):
        outputs["kleborate"] = "results/amr_detect/{sample}/kleborate.tsv"
    elif "Staphylococcus" in taxa and "aureus" in taxa:
        outputs["spa_typer"] = "results/amr_detect/{sample}/spa_typer.tsv"
        outputs["SCCmec"] = "results/amr_detect/{sample}/SCCmec.tsv"
    elif taxa.startswith("Escherichia") or taxa.startswith("Shigella"):
        outputs["etoki_ebeis"] = "results/amr_detect/{sample}/etoki_ebeis.tsv"
    elif taxa.startswith("Salmonella"):
        outputs["sistr"] = "results/amr_detect/{sample}/sistr_serovar.tab"
        outputs["seroseq"] = "results/amr_detect/{sample}/seqsero_summary.tsv"
        if "enterica" in taxa:
            outputs["crispol"] = "results/amr_detect/{sample}/crispol.tsv"
    elif taxa.startswith("Streptococcus") and "pneumoniae" in taxa:
        outputs["pneumokity"] = "results/amr_detect/{sample}/pneumo_capsular_typing/{sample}_result_data.csv"
    try:
        matched_organism = get_key_for_value_from_db(taxa, MLST_MAP)
        outputs["mlst"] = "results/amr_detect/{sample}/mlst.tsv"
        if find_cc_profile_for_taxonomy(taxa):
            outputs["clonal_complex"] = "results/amr_detect/{sample}/clonal_complex.tsv"
        else:
            logger.warning(f"Skipping clonal complex profiling for {taxa =} and {sample =} as no profile found")
    except KeyError:
        pass

    return outputs


def get_taxonomy_dependant_custom_mlst_scheme(taxa: str) -> str:
    if taxa.startswith("Bordetella") and "pertussis" in taxa:
        return "bpertussis_antigens"
    return None


def infer_outputs_for_sample(wildcards) -> dict[str, str]:
    if check_all_checks_success_for_sample(wildcards.sample):
        taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)

        outputs = {
            "taxonomy": "results/taxonomy/{sample}/parsed_taxa.txt",
            "ncbi_taxonomy_id": "results/taxonomy/{sample}/ncbi_taxa.tsv",
            "amrfinder": "results/amr_detect/{sample}/amrfinder.tsv",
            "abricate": "results/amr_detect/{sample}/abricate.tsv",
            "rgi": "results/amr_detect/{sample}/rgi_main.txt",
            "resfinder": "results/amr_detect/{sample}/resfinder/ResFinder_results_tab.txt",
            "pointfinder": "results/amr_detect/{sample}/resfinder/PointFinder_results.txt",
            "plasmids": "results/plasmids/{sample}/mob_typer.txt",
            "qc_checks": "results/checks/{sample}/qc_summary.tsv",
            "others_ncbi_tax": "results/taxonomy/{sample}/ncbi_gtdb_mapping.tsv",
        }
        if config["run_hamronization"]:
            outputs["hamronization"] = "results/hamronization/summary/{sample}.tsv"

        if not sample_has_asssembly_as_input(wildcards.sample):
            outputs["seqkit"] = "results/assembly/{sample}/seqkit_stats.tsv"

        custom_mlst_scheme = get_taxonomy_dependant_custom_mlst_scheme(taxa)
        if custom_mlst_scheme is not None:
            outputs["mlst_custom"] = "results/amr_detect/{sample}/mlst_custom.tsv"

        taxa_outputs = get_taxonomy_dependant_outputs(wildcards.sample, taxa)
        return outputs | taxa_outputs
    else:
        return {
            "qc_checks": "results/checks/{sample}/qc_summary.tsv",
        }


def infer_outputs_for_sample_as_list(wildcards):
    return list(infer_outputs_for_sample(wildcards).values())


def infer_results_to_summarize_for_sample(wildcards):
    dct = infer_outputs_for_sample(wildcards)
    reports = [
        "taxonomy",
        "ncbi_taxonomy_id",
        "mlst",
        "mlst_custom",
        "clonal_complex",
        "spa_typer",
        "SCCmec",
        "sistr",
        "seroseq",
        "kleborate",
        "plasmids",
        "amrfinder",
        "abricate",
        "seqkit",
        "qc_checks",
        "pneumokity",
        "etoki_ebeis",
    ]
    if sample_has_asssembly_as_input(wildcards.sample):
        reports.remove("seqkit")

    out_dict = {}
    for report in reports:
        if report in dct:
            out_dict[report] = dct[report]
    return out_dict


def get_all_clonal_complex_profiles():
    profiles = []
    for profile in list(config["clonal_complex"]["mapping_to_gtdbtk_names"].keys()):
        profiles.append(os.path.join(config["clonal_complex"]["db_dir"], f"{profile}.tsv"))
    return profiles


### Wildcard handling #################################################################################################


def infer_assembly_fasta(wildcards) -> str:
    if sample_has_asssembly_as_input(wildcards.sample):
        return get_fasta_for_sample_from_pep(wildcards.sample)
    else:
        return "results/assembly/{sample}/assembly_cleaned.fasta"


def infer_reads_for_assembly(wildcards) -> dict[str, str]:
    inputs = {}
    if sample_has_long_reads(wildcards.sample):
        if config["assembly__unicycler"]["use_long_if_relevant"]:
            inputs["long"] = get_long_reads_for_sample_from_pep(wildcards.sample)
        else:
            logger.warning(
                f"Long reads are available for sample={wildcards.sample} but they are ignored for assembly, as config assembly__unicycler->use_long_if_relevant is False",
            )
    inputs["paired"] = ["results/reads/trimmed/{sample}_R1.fastq.gz", "results/reads/trimmed/{sample}_R2.fastq.gz"]
    return inputs


def infer_fastqs_for_trimming(wildcards) -> list[str]:
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


def infer_resfinder_input(wildcards):
    if config["resfinder"]["input_to_use"] == "assembly":
        return infer_assembly_fasta(wildcards)
    if sample_has_asssembly_as_input(wildcards.sample):
        raise ValueError(
            f"resfinder input_to_use is set to reads, but pepfile gives that input for sample={wildcards.sample} is assembly not reads"
        )
    return ["results/reads/trimmed/{sample}_R1.fastq.gz", "results/reads/trimmed/{sample}_R2.fastq.gz"]


def infer_fastq_path_for_fastqc(wildcards):
    if wildcards.step != "original":
        return "results/reads/{step}/{sample}_{pair}.fastq.gz"
    if "pair" not in wildcards or wildcards.pair == "R1":
        return get_first_fastq_for_sample_from_pep(wildcards.sample, read_pair="fq1")
    elif wildcards.pair == "R2":
        return get_first_fastq_for_sample_from_pep(wildcards.sample, read_pair="fq2")


def get_organism_for_amrfinder(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)
    try:
        matched_organism = get_key_for_value_from_db(taxa, AMRFINDER_MAP)
        return f"--organism {matched_organism}"
    except KeyError:
        logger.warning(f"Could not find organism {taxa} for sample {wildcards.sample} in amrfinder map")
        return ""


def get_taxonomy_for_mlst(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)
    try:
        matched_organism = get_key_for_value_from_db(taxa, MLST_MAP)
        return f"--scheme {matched_organism}"
    except KeyError:
        logger.warning(f"Could not find organism {taxa} for sample {wildcards.sample} in MLST map")
        return ""


def get_custom_scheme_for_mlst(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)
    scheme = get_taxonomy_dependant_custom_mlst_scheme(taxa)
    return f"--scheme {scheme}"


def find_cc_profile_for_taxonomy(taxa: str):
    for profile, organism_regex in config["clonal_complex"]["mapping_to_gtdbtk_names"].items():
        if re.match(organism_regex, taxa):
            return profile
    return None


def infer_profile_for_clonal_complex(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)
    profile = find_cc_profile_for_taxonomy(taxa)
    return os.path.join(config["clonal_complex"]["db_dir"], f"{profile}.tsv")


def get_taxonomy_for_resfinder(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)
    if len(taxa) < 2:
        return "unknown species"
    return taxa.lower()


def infer_relevant_checks(wildcards):
    if sample_has_asssembly_as_input(wildcards.sample):
        return ["results/checks/{sample}/check_skipping.tsv"]

    checks = ["results/checks/{sample}/pre_assembly_summary.tsv"]

    if not check_preassembly_QC_for_sample(wildcards.sample):
        return checks

    checks += [
        "results/checks/{sample}/assembly_constructed.tsv",
    ]

    if check_assembly_construction_success_for_sample(wildcards.sample):
        checks += [
            "results/checks/{sample}/assembly_quality.tsv",
            "results/checks/{sample}/coverage_check.tsv",
            "results/checks/{sample}/self_contamination_check.tsv",
        ]

    return checks


### Parameter parsing #################################################################################################


def get_unicycler_params():
    extra = f"--min_fasta_length {config['assembly__unicycler']['min_fasta_length']}"
    extra += f" --mode {config['assembly__unicycler']['bridging_mode']}"
    extra += f" --linear_seqs {config['assembly__unicycler']['linear_seqs']}"
    return extra


def parse_adapter_removal_params(adapter_config):
    args_lst = [
        f"--action {adapter_config['action']}",
        f"--overlap {adapter_config['overlap']}",
        f"--times {adapter_config['times']}",
        f"--error-rate {adapter_config['error_rate']}",
    ]
    if adapter_config["keep_trimmed_only"]:
        args_lst.append("--discard-untrimmed")

    if path := adapter_config["adapters_anywhere_file"]:
        if not os.path.exists(path):
            raise ValueError(
                f"Adapter removal is enabled, but the adapter file specified in the adapters_anywhere_file element = {path} does not exist."
            )
        args_lst.append(f"--anywhere file:{path} -B file:{path}")

    if end3_file := adapter_config["adapters_3_end_file"]:
        if not os.path.exists(end3_file):
            raise ValueError(
                f"Adapter removal is enabled, but the adapter file specified in the adapters_3_end_file element = {end3_file} does not exist."
            )
        args_lst.append(f"--adapter file:{path} -A file:{path}")

    if end5_file := adapter_config["adapters_5_end_file"]:
        if not os.path.exists(end5_file):
            raise ValueError(
                f"Adapter removal is enabled, but the adapter file specified in the adapters_5_end_file element = {end5_file} does not exist."
            )
        args_lst.append(f"--front file:{path} -G file:{path}")

    return args_lst


def get_cutadapt_extra(cutadapt_config) -> list[str]:
    args_lst = []

    if (value := cutadapt_config["shorten_to_length"]) is not None:
        args_lst.append(f"--length {value}")
    if (value := cutadapt_config["cut_from_start_r1"]) is not None:
        args_lst.append(f"--cut {value}")
    if (value := cutadapt_config["cut_from_start_r2"]) is not None:
        args_lst.append(f"-U {value}")
    if (value := cutadapt_config["cut_from_end_r1"]) is not None:
        args_lst.append(f"--cut -{value}")
    if (value := cutadapt_config["cut_from_end_r2"]) is not None:
        args_lst.append(f"-U -{value}")

    if (value := cutadapt_config["max_n_bases"]) is not None:
        args_lst.append(f"--max-n {value}")
    if (value := cutadapt_config["max_expected_errors"]) is not None:
        args_lst.append(f"--max-expected-errors {value}")
    if value := cutadapt_config["trim_N_bases_on_ends"]:
        args_lst.append(f"--trim-n")
    if cutadapt_config["nextseq_trimming_mode"]:
        value = config["reads__trimming"]["quality_cutoff_from_3_end_r1"]
        args_lst.append(f"--nextseq-trim={value}")

    if cutadapt_config["do_adapter_removal"]:
        args_lst += parse_adapter_removal_params(cutadapt_config["adapter_removal"])

    return args_lst


def parse_paired_cutadapt_param(pe_config, param1, param2, arg_name) -> str:
    if pe_config.get(param1, None) is not None:
        if pe_config.get(param2, None) is not None:
            return f"{arg_name} {pe_config[param1]}:{pe_config[param2]}"
        else:
            return f"{arg_name} {pe_config[param1]}:"
    elif pe_config.get(param2, None) is not None:
        return f"{arg_name} :{pe_config[param2]}"
    return ""


def parse_cutadapt_comma_param(cutadapt_config, param1, param2, arg_name) -> str:
    if cutadapt_config.get(param1) is not None:
        if cutadapt_config.get(param2) is not None:
            return f"{arg_name} {cutadapt_config[param2]},{cutadapt_config[param1]}"
        else:
            return f"{arg_name} {cutadapt_config[param1]}"
    elif cutadapt_config.get(param2) is not None:
        return f"{arg_name} {cutadapt_config[param2]},0"
    return ""


def get_cutadapt_extra_pe() -> str:
    cutadapt_config = config["reads__trimming"]

    args_lst = get_cutadapt_extra(cutadapt_config)

    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "max_length_r1", "max_length_r2", "--maximum-length"):
        args_lst.append(parsed_arg)
    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "min_length_r1", "min_length_r2", "--minimum-length"):
        args_lst.append(parsed_arg)
    if qual_cut_arg_r1 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "--quality-cutoff"
    ):
        args_lst.append(qual_cut_arg_r1)
    if qual_cut_arg_r2 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "-Q"
    ):
        args_lst.append(qual_cut_arg_r2)
    return " ".join(args_lst)


### Resource handling #################################################################################################


def get_mem_mb_for_trimming(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["reads__trimming_mem_mb"] * attempt)


def get_mem_mb_for_fastqc(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["fastqc_mem_mb"] * attempt)


def get_mem_mb_for_unicycler(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["assembly__unicycler_mem_mb"] * attempt)


def get_mem_mb_for_gtdb(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["gtdb_classify__mem_mb"] * attempt)


def get_mem_mb_for_mapping_postprocess(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__mem_mb"] * attempt)


def get_mem_mb_for_mapping(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping_postprocess__mem_mb"] * attempt)
