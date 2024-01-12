from snakemake.utils import validate


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
        args_lst.append(f"--cut {config['reads__trimming']['cut_from_start']}")
    if "cut_from_end" in config["reads__trimming"]:
        args_lst.append(f"--cut -{config['reads__trimming']['cut_from_end']}")
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


def optional_bandage_outputs(wildcards):
    if check_assembly_construction_success_for_sample(wildcards.sample):
        return [
            "results/assembly/{sample}/bandage/bandage.svg",
            "results/assembly/{sample}/bandage/bandage.info",
            "results/taxonomy/{sample}",
        ]
    else:
        return "results/checks/{sample}/assembly_constructed.txt"


def bandage_check_if_relevant(wildcards):
    if check_assembly_construction_success_for_sample(wildcards.sample):
        return "results/checks/{sample}/assembly_quality.txt"
    else:
        return ""


def check_assembly_construction_success_for_sample(sample: str):
    with checkpoints.assembly_constructed.get(sample=sample).output[0].open() as f:
        return f.read().startswith("PASS:")


def check_all_checks_success_for_sample(sample: str):
    with checkpoints.summary_all_checks.get(sample=sample).output[0].open() as f:
        return all([line.startswith("PASS:") for line in f.readlines()])


def get_second_phase_results(wildcards):
    base_result = ["results/checks/{sample}/summary.txt"]

    if check_all_checks_success_for_sample(wildcards.sample):
        base_result.append("results/amr_detect/{sample}/amrfinder")
        base_result.append("results/amr_detect/{sample}/mlst.tsv")

    return base_result


def get_all_checks(wildcards):
    basic_checks = [
        "results/checks/{sample}/foreign_contamination.txt",
        "results/checks/{sample}/assembly_constructed.txt",
    ]

    if check_assembly_construction_success_for_sample(wildcards.sample):
        basic_checks.append("results/checks/{sample}/assembly_quality.txt")

    return basic_checks


def get_outputs():
    sample_names = get_sample_names()
    return {
        "bandage_reports_optional": expand("results/checks/{sample}/.bandage_requested.txt", sample=sample_names),
        "multiqc": expand("results/multiqc/{sample}.html", sample=sample_names),
        "checks": expand("results/checks/{sample}/.final_results_requested.txt", sample=sample_names),
    }


### Resource handling #################################################################################################


def get_mem_mb_for_trimming(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["reads__trimming_mem_mb"] * attempt)


def get_mem_mb_for_fastqc(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["fastqc_mem_mb"] * attempt)


def get_mem_mb_for_unicycler(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["assembly__unicycler_mem_mb"] * attempt)


def get_mem_mb_for_gtdb(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["gtdb_classify__mem_mb"] * attempt)
