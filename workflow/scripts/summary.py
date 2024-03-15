import functools
import os
import sys


def parse_amrfinder(path: str) -> dict[str, str]:
    outputs: dict[tuple[str, str], str] = {}
    tmp_outputs: dict[tuple[str, str], list[str]] = {}
    with open(path, "r") as f:
        header = f.readline().rstrip().split("\t")
        gene_symbol_idx = header.index("Gene symbol")
        class_idx = header.index("Class")
        subclass_idx = header.index("Subclass")
        rows = [row.rstrip().split("\t") for row in f.readlines()]
        for row in rows:
            if (row[class_idx], row[subclass_idx]) in tmp_outputs:
                tmp_outputs[(row[class_idx], row[subclass_idx])].append(row[gene_symbol_idx])
            else:
                tmp_outputs[(row[class_idx], row[subclass_idx])] = [row[gene_symbol_idx]]

        for k, v in tmp_outputs.items():
            sorted_vals = sorted(v)
            outputs[k] = ", ".join(sorted_vals)

    sorted(outputs)
    return {":AMRFINDER:".join(x): y for x, y in outputs.items()}


def index_based_parser(path: str, indexes: list[int], recode_into_columns: list[str]):
    with open(path, "r") as f:
        row = f.readline().rstrip().split("\t")
        return {col: row[idx] for col, idx in zip(recode_into_columns, indexes)}


def column_based_parser(path: str, columns: list[str], recode_into_columns: list[str] | None = None):
    if not recode_into_columns:
        recode_into_columns = columns
    with open(path, "r") as f:
        header = f.readline().rstrip().split("\t")
        row = f.readline().rstrip().split("\t")
        return {new_col: row[header.index(col)] for col, new_col in zip(columns, recode_into_columns)}


def transpose_multi_columns_parser(path: str, transpose: tuple[str, str]):
    out: dict[str, str] = {}
    with open(path, "r") as f:
        header = f.readline().rstrip().split("\t")
        transposed_idx = header.index(transpose[0])
        value_idx = header.index(transpose[1])
        rows = [row.rstrip().split("\t") for row in f.readlines()]
        if len(rows) == 1 and "assembly_not_requested" not in rows[0]:
            out["running_from_assembly"] = "PASS"
            return out
        for row in rows:
            out[row[transposed_idx]] = row[value_idx]
        return out


def row_joiner_on_column_parser(
    path: str, join_multiple_rows_on: str, columns: list[str], recode_into_columns: list[str] | None = None
):
    if not recode_into_columns:
        recode_into_columns = columns
    with open(path, "r") as f:
        header = f.readline().rstrip().split("\t")
        rows = [row.rstrip().split("\t") for row in f.readlines()]
        out: dict[str, str] = {}
        for column, new_column in zip(columns, recode_into_columns):
            idx = header.index(column)
            vals = join_multiple_rows_on.join([row[idx] for row in rows])
            out[new_column] = vals
        return out


mapping_functions = {
    "taxonomy": functools.partial(index_based_parser, indexes=[0], recode_into_columns=["taxonomy"]),
    "mlst": functools.partial(index_based_parser, indexes=[2], recode_into_columns=["mlst"]),
    "spatyper": functools.partial(column_based_parser, columns=["Type"], recode_into_columns=["spa_type"]),
    "sistr": functools.partial(column_based_parser, columns=["h1", "h2", "o_antigen", "serogroup", "serovar"]),
    "seroseq": functools.partial(column_based_parser, columns=["Predicted antigenic profile", "Predicted serotype"]),
    "qc_checks": functools.partial(transpose_multi_columns_parser, transpose=("parameter", "result")),
    "kleborate": functools.partial(
        column_based_parser,
        columns=[
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
        ],
    ),
    "plasmids": functools.partial(column_based_parser, columns=["rep_type(s)", "predicted_mobility"]),
    "seqkit": functools.partial(column_based_parser, columns=["num_seqs", "sum_len"]),
    "amrfinder": parse_amrfinder,
    "abricate": functools.partial(
        row_joiner_on_column_parser,
        join_multiple_rows_on=", ",
        columns=["GENE"],
        recode_into_columns=["virulence_database_hits"],
    ),
}


def run(
    results: dict[str, str], output_file: str, out_delimiter: str, missing_values_placeholder: str, sample_name: str
):

    print("Processing results...", file=sys.stderr)
    print(f"Results: {results}", file=sys.stderr)
    header_value_dict: dict[str, str] = {"sample": sample_name}
    for result_type, result_path in results.items():
        header_value_dict = header_value_dict | mapping_functions[result_type](result_path)

    header = out_delimiter.join(header_value_dict.keys())
    values = out_delimiter.join(header_value_dict.values())

    if not os.path.isdir(os.path.dirname(output_file)):
        print("Creating output directory...", file=sys.stderr)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        f.write(f"{header}\n{values}\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    run(
        results=snakemake.input,
        output_file=snakemake.output.tsv,
        out_delimiter=snakemake.params.delimiter,
        missing_values_placeholder=snakemake.params.placeholder,
        sample_name=snakemake.params.sample_name,
    )
