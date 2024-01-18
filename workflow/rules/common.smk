from snakemake.utils import validate
import json
import re


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")


pepfile: config["pepfile"]


validate(pep.sample_table, "../schemas/samples.schema.yaml")


def get_sample_names():
    return list(pep.sample_table["sample_name"].values)


def get_constraints():
    constraints = {
        "sample": "|".join(get_sample_names()),
    }
    return constraints


def get_reads_for_trimming(wildcards):
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


def get_one_fastq_file(sample: str, read_pair="fq1"):
    return pep.sample_table.loc[sample][[read_pair]]


def infer_fastq_path(wildcards):
    if wildcards.step != "original":
        return "results/reads/{step}/{sample}_{pair}.fastq.gz"
    if "pair" not in wildcards or wildcards.pair == "R1":
        return get_one_fastq_file(wildcards.sample, read_pair="fq1")[0]
    elif wildcards.pair == "R2":
        return get_one_fastq_file(wildcards.sample, read_pair="fq2")[0]


with open(f"{workflow.basedir}/resources/gtdb_amrfinder.json", "r") as f:
    AMRFINDER_MAP = json.load(f)
with open(f"{workflow.basedir}/resources/gtdb_mlst.json", "r") as f:
    MLST_MAP = json.load(f)

### Data input handling independent of wildcards ######################################################################


def get_unicycler_params():
    extra = f"--min_fasta_length {config['assembly__unicycler']['min_fasta_length']}"
    extra += f" --mode {config['assembly__unicycler']['bridging_mode']}"
    extra += f" --linear_seqs {config['assembly__unicycler']['linear_seqs']}"
    return extra


def parse_adapter_removal_params():
    args_lst = []
    adapters_file = config["reads__trimming"]["adapter_removal"]["adapters_fasta"]
    read_location = config["reads__trimming"]["adapter_removal"]["read_location"]
    args_lst.append(f"--{read_location} file:{adapters_file}")

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
    if "cut_from_start" in config["reads__trimming"]:
        args_lst.append(
            f"--cut {config['reads__trimming']['cut_from_start']} -U {config['reads__trimming']['cut_from_start']}"
        )
    if "cut_from_end" in config["reads__trimming"]:
        args_lst.append(
            f"--cut -{config['reads__trimming']['cut_from_end']} -U -{config['reads__trimming']['cut_from_end']}"
        )
    if "max_n_bases" in config["reads__trimming"]:
        args_lst.append(f"--max-n {config['reads__trimming']['max_n_bases']}")
    if "max_expected_errors" in config["reads__trimming"]:
        args_lst.append(f"--max-expected-errors {config['reads__trimming']['max_expected_errors']}")

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


### Global rule-set stuff #############################################################################################


def post_assembly_outputs(wildcards):
    if check_assembly_construction_success_for_sample(wildcards.sample) and not config["gtdb_hack"]:
        return [
            "results/assembly/{sample}/bandage/bandage.svg",
            "results/assembly/{sample}/bandage/bandage.info",
            "results/taxonomy/{sample}",
        ]
    else:
        return "results/checks/{sample}/assembly_constructed.txt"


def check_assembly_construction_success_for_sample(sample: str):
    with checkpoints.assembly_constructed.get(sample=sample).output[0].open() as f:
        return f.read().startswith("PASS:")


def check_all_checks_success_for_sample(sample: str):
    with checkpoints.summary_all_checks.get(sample=sample).output[0].open() as f:
        return all([line.startswith("PASS:") for line in f.readlines()])


def get_all_checks(wildcards):
    basic_checks = [
        "results/checks/{sample}/foreign_contamination.txt",
        "results/checks/{sample}/assembly_constructed.txt",
    ]

    if check_assembly_construction_success_for_sample(wildcards.sample) and not config["gtdb_hack"]:
        basic_checks.append("results/checks/{sample}/assembly_quality.txt")
        basic_checks.append("results/checks/{sample}/coverage_check.txt")
        basic_checks.append("results/checks/{sample}/self_contamination_check.txt")

    return basic_checks


def get_outputs():
    sample_names = get_sample_names()
    return {
        "bandage_reports_optional": expand("results/checks/{sample}/.bandage_requested.txt", sample=sample_names),
        "multiqc": "results/summary/multiqc.html",
        "checks": expand("results/checks/{sample}/.final_results_requested.txt", sample=sample_names),
    }


def get_parsed_taxa_from_gtdbtk_for_sample(sample: str):
    with checkpoints.gtdbtk__parse_taxa.get(sample=sample).output[0].open() as f:
        taxa = f.read().strip()
    return taxa


def get_key_for_value_from_db(value: str, db: dict):
    for key in db:
        pattern = "bob"
        if re.match(pattern, taxa):
            return key
    raise KeyError


def get_organism_for_amrfinder(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk(wildcards.sample)
    try:
        return get_key_for_value_from_db(taxa, AMRFINDER_MAP)
    except KeyError:
        raise KeyError(f"Could not find organism {taxa} for sample {wildcards.sample} in amrfinder map")


def get_taxonomy_for_mlst(wildcards):
    taxa = get_parsed_taxa_from_gtdbtk(wildcards.sample)
    try:
        return get_key_for_value_from_db(taxa, MLST_MAP)
    except KeyError:
        raise KeyError(f"Could not find organism {taxa} for sample {wildcards.sample} in MLST map")


def get_second_phase_results(wildcards):
    base_result = ["results/checks/{sample}/summary.txt"]

    if not config["gtdb_hack"] and check_all_checks_success_for_sample(wildcards.sample):
        base_result.append("results/amr_detect/{sample}/amrfinder.tsv")
        base_result.append("results/amr_detect/{sample}/mlst.tsv")
        base_result.append("results/amr_detect/{sample}/abricate.tsv")

        taxa = get_parsed_taxa_from_gtdbtk_for_sample(wildcards.sample)
        if taxa.startswith("Klebsiella"):
            base_result.append("results/amr_detect/{sample}/kleborate.tsv")
        elif "Staphylococcus" in taxa and "aureus" in taxa:
            base_result.append("results/amr_detect/{sample}/spa_typer.tsv")
        elif taxa.startswith("Escherichia") or taxa.startswith("Shigella"):
            base_result.append("results/amr_detect/{sample}/etoki_ebeis.tsv")
    return base_result


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
