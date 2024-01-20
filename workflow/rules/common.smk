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
        assembly = pep.sample_table.loc[sample][["fasta"]]
        return True
    except KeyError:
        return False


def get_sample_names_with_reads_as_input() -> list[str]:
    return [sample for sample in get_sample_names() if not sample_has_asssembly_as_input(sample)]


def get_constraints() -> dict[str, list[str]]:
    return {"sample": "|".join(get_sample_names()), "pair": ["R1", "R2"]}


def get_first_fastq_for_sample_from_pep(sample: str, read_pair="fq1") -> str:
    return pep.sample_table.loc[sample][[read_pair]][0]


def get_fasta_for_sample_from_pep(sample: str) -> str:
    return pep.sample_table.loc[sample][["fasta"]][0]


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
        if re.match(db[key]["genus"], value):
            return key
    raise KeyError


def get_parsed_taxa_from_gtdbtk_for_sample(sample: str):
    with checkpoints.gtdbtk__parse_taxa.get(sample=sample).output[0].open() as f:
        return f.read().strip()


def check_assembly_construction_success_for_sample(sample: str):
    with checkpoints.assembly_constructed.get(sample=sample).output[0].open() as f:
        return f.read().startswith("PASS:")


def check_all_checks_success_for_sample(sample: str):
    with checkpoints.summary_all_checks.get(sample=sample).output[0].open() as f:
        return all([line.startswith("PASS:") for line in f.readlines()])


def get_outputs():
    sample_names = get_sample_names()
    outputs = {
        "final_results": expand("results/checks/{sample}/.final_results_requested.txt", sample=sample_names),
    }

    if len(get_sample_names_with_reads_as_input()) > 1:
        outputs["multiqc"] = "results/summary/multiqc.html"

    return outputs


def get_taxonomy_dependant_outputs(sample: str, taxa: str) -> list[str]:
    outputs = []
    if taxa.startswith("Klebsiella"):
        outputs.append("results/amr_detect/{sample}/kleborate.tsv")
    elif "Staphylococcus" in taxa and "aureus" in taxa:
        outputs.append("results/amr_detect/{sample}/spa_typer.tsv")
    elif taxa.startswith("Escherichia") or taxa.startswith("Shigella"):
        outputs.append("results/amr_detect/{sample}/etoki_ebeis.tsv")
    return outputs


def infer_outputs_for_sample(wildcards):
    if check_all_checks_success_for_sample(wildcards.sample):
        taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)

        return [
            "results/amr_detect/{sample}/amrfinder.tsv",
            "results/amr_detect/{sample}/mlst.tsv",
            "results/amr_detect/{sample}/abricate.tsv",
        ] + get_taxonomy_dependant_outputs(wildcards.sample, taxa)

    else:
        return "results/checks/{sample}/summary.txt"


### Wildcard handling #################################################################################################


def infer_assembly_fasta(wildcards) -> str:
    if sample_has_asssembly_as_input(wildcards.sample):
        return get_fasta_for_sample_from_pep(wildcards.sample)
    else:
        return "results/assembly/{sample}/assembly.fasta"


def infer_fastqs_for_trimming(wildcards) -> list[str]:
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


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
        return get_key_for_value_from_db(taxa, AMRFINDER_MAP)
    except KeyError:
        raise KeyError(f"Could not find organism {taxa} for sample {wildcards.sample} in amrfinder map")


def get_taxonomy_for_mlst(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)
    try:
        return get_key_for_value_from_db(taxa, MLST_MAP)
    except KeyError:
        raise KeyError(f"Could not find organism {taxa} for sample {wildcards.sample} in MLST map")


def infer_relevant_checks(wildcards):
    if sample_has_asssembly_as_input(wildcards.sample):
        return ["results/checks/{sample}/check_skipping.txt"]

    checks = [
        "results/checks/{sample}/foreign_contamination.txt",
        "results/checks/{sample}/assembly_constructed.txt",
    ]

    if check_assembly_construction_success_for_sample(wildcards.sample) and not config["gtdb_hack"]:
        checks += [
            "results/checks/{sample}/assembly_quality.txt",
            "results/checks/{sample}/coverage_check.txt",
            "results/checks/{sample}/self_contamination_check.txt",
        ]

    return checks


### Parameter parsing #################################################################################################


def get_unicycler_params():
    extra = f"--min_fasta_length {config['assembly__unicycler']['min_fasta_length']}"
    extra += f" --mode {config['assembly__unicycler']['bridging_mode']}"
    extra += f" --linear_seqs {config['assembly__unicycler']['linear_seqs']}"
    return extra


def parse_adapter_removal_params():
    args_lst = []
    adapters_file = config["reads__trimming"]["adapter_removal"]["adapters_fasta"]
    read_location = config["reads__trimming"]["adapter_removal"]["read_location"]

    if read_location == "front":
        paired_arg = "-G"
    elif read_location == "anywhere":
        paired_arg = "-B"
    elif read_location == "adapter":
        paired_arg = "-A"

    args_lst.append(f"--{read_location} file:{adapters_file} {paired_arg} file:{adapters_file}")

    if config["reads__trimming"]["adapter_removal"]["keep_trimmed_only"]:
        args_lst.append("--discard-untrimmed")

    args_lst.append(f"--action {config['reads__trimming']['adapter_removal']['action']}")
    args_lst.append(f"--overlap {config['reads__trimming']['adapter_removal']['overlap']}")
    args_lst.append(f"--times {config['reads__trimming']['adapter_removal']['times']}")
    args_lst.append(f"--error-rate {config['reads__trimming']['adapter_removal']['error_rate']}")
    return args_lst


def get_cutadapt_extra() -> list[str]:
    args_lst = []

    if "shorten_to_length" in config["reads__trimming"]:
        args_lst.append(f"--length {config['reads__trimming']['shorten_to_length']}")
    if value := config["reads__trimming"].get("cut_from_start_r1", None):
        args_lst.append(f"--cut {value}")
    if value := config["reads__trimming"].get("cut_from_start_r2", None):
        args_lst.append(f"-U {value}")
    if value := config["reads__trimming"].get("cut_from_end_r1", None):
        args_lst.append(f"--cut -{value}")
    if value := config["reads__trimming"].get("cut_from_end_r2", None):
        args_lst.append(f"-U -{value}")

    if "max_n_bases" in config["reads__trimming"]:
        args_lst.append(f"--max-n {config['reads__trimming']['max_n_bases']}")
    if "max_expected_errors" in config["reads__trimming"]:
        args_lst.append(f"--max-expected-errors {config['reads__trimming']['max_expected_errors']}")
    if "trim_N_bases_on_ends" in config["reads__trimming"]:
        args_lst.append(f"--trim-n")

    if config["reads__trimming"]["adapter_removal"]["do"]:
        args_lst += parse_adapter_removal_params()

    return args_lst


def parse_paired_cutadapt_param(pe_config, param1, param2, arg_name) -> str:
    if param1 in pe_config:
        if param2 in pe_config:
            return f"{arg_name} {pe_config[param1]}:{pe_config[param2]}"
        else:
            return f"{arg_name} {pe_config[param1]}:"
    elif param2 in pe_config:
        return f"{arg_name} :{pe_config[param2]}"
    return ""


def parse_cutadapt_comma_param(config, param1, param2, arg_name) -> str:
    if param1 in config:
        if param2 in config:
            return f"{arg_name} {config[param2]},{config[param1]}"
        else:
            return f"{arg_name} {config[param1]}"
    elif param2 in config:
        return f"{arg_name} {config[param2]},0"
    return ""


def get_cutadapt_extra_pe() -> str:
    args_lst = get_cutadapt_extra()

    cutadapt_config = config["reads__trimming"]
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
