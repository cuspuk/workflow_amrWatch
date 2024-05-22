import os
import sys


def load_tsv(tsv_path: str) -> dict[str, str]:
    with open(tsv_path, "r") as f:
        header = f.readline().strip().split("\t")
        row = f.readline().strip().split("\t")
        if len(header) != len(row):
            raise ValueError("Header and row length mismatch")
        return dict(zip(header, row))


def run(tsvs: list[str], output_file: str, out_delimiter: str, nan_value: str, amrfinder_uniq_tag: str):

    analysis_results: list[dict[str, str]] = [load_tsv(tsv) for tsv in tsvs]

    base_columns = [
        "sample",
        "taxonomy",
        "ncbi_taxonomy_id",
        "mlst",
        "serotype",
        "clonal_complex",
        "spa_type",
        "SCCmec_type",
        "SCCmecA_presence",
    ]
    qc_columns = [
        "assembly_length",
        "number_of_contigs",
        "number_of_dead_ends",
        "mean_coverage",
        "foreign_contamination",
        "trimmed_fastq_length",
        "ambiguous_bases",
        "human_contamination",
        "assembly_construction",
        "assembly_not_requested",
    ]
    full_qc_columns: list[str] = []
    for col in qc_columns:
        full_qc_columns.extend([f"{col}{suffix}" for suffix in ["__result", "__value", "__comment"]])

    technical_columns = [
        "num_seqs",
        "sum_len",
    ]
    salmonella_columns = [
        "h1",
        "h2",
        "o_antigen",
        "serogroup",
        "serovar",
        "Predicted antigenic profile",
        "Predicted serotype",
    ]
    kleborate_columns = [
        "virulence_score",
        "resistance_score",
        "num_resistance_classes",
        "num_resistance_genes",
        "Yersiniabactin",
        "YbST",
        "Colibactin",
        "CbST",
        "Aerobactin",
        "AbST",
        "Salmochelin",
        "SmST",
        "RmpADC",
        "RmST",
        "rmpA2",
        "wzi",
    ]
    plasmids_columns = ["rep_type(s)", "predicted_mobility"]
    abricate_columns = ["virulence_database_hits"]

    amrfinder_tuples: list[tuple[str, str]] = []
    for result in analysis_results:
        tuples = [
            (x.split(amrfinder_uniq_tag)[0], x.split(amrfinder_uniq_tag)[1])
            for x in result.keys()
            if amrfinder_uniq_tag in x and "__result" not in x and "__value" not in x and "__comment" not in x
        ]
        amrfinder_tuples.extend(
            [amrfinder_tuple for amrfinder_tuple in tuples if amrfinder_tuple not in amrfinder_tuples]
        )

    amrfinder_tuples.sort()
    amrfinder_columns = [f"{x}{amrfinder_uniq_tag}{y}" for x, y in amrfinder_tuples]

    common_header: list[str] = (
        base_columns
        + full_qc_columns
        + technical_columns
        + salmonella_columns
        + kleborate_columns
        + plasmids_columns
        + amrfinder_columns
        + abricate_columns
    )

    merged_results: list[dict[str, str]] = [
        {column: result.get(column, nan_value) for column in common_header} for result in analysis_results
    ]
    missing_keys = [
        x
        for result in analysis_results
        for x in result.keys()
        if x not in common_header and amrfinder_uniq_tag not in x
    ]
    if missing_keys:
        raise ValueError(f"There have been some keys that were not contained in summary: {missing_keys=}")

    if not os.path.isdir(os.path.dirname(output_file)):
        print("Creating output directory...", file=sys.stderr)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        f.write(out_delimiter.join(common_header) + "\n")
        for row in merged_results:
            f.write(out_delimiter.join(list(row.values())) + "\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    run(
        tsvs=snakemake.input.tsvs,
        output_file=snakemake.output.tsv,
        out_delimiter=snakemake.params.delimiter,
        amrfinder_uniq_tag=snakemake.params.amrfinder_uniq_tag,
        nan_value=snakemake.params.nan_value,
    )
