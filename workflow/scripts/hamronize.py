from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with open(snakemake.input.version, "r") as f:
    version = f.read().strip()

with open(snakemake.input.db_version, "r") as f:
    db_version = f.read().strip()

input_file_name = ""
if snakemake.params.tool in ["rgi", "abricate"]:
    input_file_name = f"--input_file_name {snakemake.params.sample_name}"

shell(
    "hamronize {snakemake.params.tool}"
    " {input_file_name}"
    " --analysis_software_version {version}"
    " --reference_database_version {db_version}"
    " {snakemake.input.rgi} --output {snakemake.output.tsv} --format tsv {log}"
)
