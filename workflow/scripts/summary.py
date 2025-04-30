import csv
import functools
import os
import re
import sys


def atoi(text: str):
    return int(text) if text.isdigit() else text.lower()


def human_sort(text: str):
    """
    Taken from http://nedbatchelder.com/blog/200712/human_sorting.html
    """
    return [atoi(c) for c in re.split(r"(\d+)", text)]


def parse_amrfinder(path: str, amrfinder_uniq_tag: str) -> dict[str, str]:
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
    return {amrfinder_uniq_tag.join(x): y for x, y in outputs.items()}


def index_based_parser(path: str, indexes: list[int], recode_into_columns: list[str]):
    with open(path, "r") as f:
        row = f.readline().rstrip().split("\t")
        return {col: row[idx] for col, idx in zip(recode_into_columns, indexes)}


def column_based_parser(path: str, columns: list[str], recode_into_columns: list[str] | None = None, delim: str = "\t"):
    if not recode_into_columns:
        recode_into_columns = columns
    with open(path, "r") as f:
        csvreader = csv.reader(f, delimiter=delim, quotechar='"')
        header = next(csvreader)
        row = next(csvreader)
        return {new_col: row[header.index(col)] for col, new_col in zip(columns, recode_into_columns)}


def etoki_ebeis_parser(path: str, columns: list[str], recode_into_columns: list[str] | None = None):
    if not recode_into_columns:
        recode_into_columns = columns

    import json

    with open(path) as f:
        data = json.load(f)

    return {new_col: data[column_name] for column_name, new_col in zip(columns, recode_into_columns)}


def sccmec_parser(
    path: str,
    columns: list[str],
    merged_column: str,
    ignore_columns: list[str],
    recode_into_columns: list[str] | None = None,
):
    if not recode_into_columns:
        recode_into_columns = columns
    with open(path, "r") as f:
        header = f.readline().rstrip().split("\t")
        row = f.readline().rstrip().split("\t")
    outs = {new_col: row[header.index(col)] for col, new_col in zip(columns, recode_into_columns)}

    outs[merged_column] = ", ".join(
        [col for col in header if col not in ignore_columns and row[header.index(col)] == "True"]
    )
    return outs


def transpose_multi_columns_parser(path: str, transpose: tuple[str, list[str]]):
    out: dict[str, str] = {}
    with open(path, "r") as f:
        header = f.readline().rstrip().split("\t")
        transposed_idx = header.index(transpose[0])
        rows = [row.rstrip().split("\t") for row in f.readlines()]

        for key in transpose[1]:
            key_idx = header.index(key)
            for row in rows:
                out[f"{row[transposed_idx]}__{key}"] = row[key_idx]
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
            vals = [row[idx] for row in rows]
            vals.sort(key=human_sort)
            merged_vals = join_multiple_rows_on.join([val for val in vals])
            out[new_column] = merged_vals
        return out


def aggregate_serotypes(serotypes: dict[str, str]) -> str | None:

    if serotypes.get("Predicted serotype", None):
        return serotypes["Predicted serotype"]

    if serotypes.get("pneumo_capsular_type", None):
        return serotypes["pneumo_capsular_type"]

    if "escherichia_o_antigen" in serotypes and "escherichia_h_antigen" in serotypes:
        combined = f"O={serotypes['escherichia_o_antigen']},H={serotypes['escherichia_h_antigen']}"
        return combined
    return None


def mlst_custom_parser(path: str):
    with open(path, "r") as f:
        header = f.readline().rstrip().split("\t")[3:]
        row = f.readline().rstrip().split("\t")[3:]
        return dict(zip(header, row))


def run(results: dict[str, str], output_file: str, out_delimiter: str, sample_name: str, amrfinder_uniq_tag: str):

    mapping_functions = {
        "taxonomy": functools.partial(index_based_parser, indexes=[0], recode_into_columns=["taxonomy"]),
        "ncbi_taxonomy_id": functools.partial(
            column_based_parser, columns=["Majority vote NCBI classification"], recode_into_columns=["ncbi_taxonomy_id"]
        ),
        "mlst": functools.partial(index_based_parser, indexes=[2], recode_into_columns=["mlst"]),
        "mlst_custom": functools.partial(mlst_custom_parser),
        "clonal_complex": functools.partial(column_based_parser, columns=["clonal_complex"]),
        "spa_typer": functools.partial(column_based_parser, columns=["Type"], recode_into_columns=["spa_type"]),
        "SCCmec": functools.partial(
            sccmec_parser,
            columns=["meca"],
            recode_into_columns=["SCCmecA_presence"],
            merged_column="SCCmec_type",
            ignore_columns=["sample", "meca"],
        ),
        "sistr": functools.partial(column_based_parser, columns=["h1", "h2", "o_antigen", "serogroup", "serovar"]),
        "seroseq": functools.partial(
            column_based_parser, columns=["Predicted antigenic profile", "Predicted serotype"]
        ),
        "etoki_ebeis": functools.partial(
            etoki_ebeis_parser,
            columns=["O", "H"],
            recode_into_columns=["escherichia_o_antigen", "escherichia_h_antigen"],
        ),
        "pneumokity": functools.partial(
            column_based_parser, columns=["predicted_serotype"], recode_into_columns=["pneumo_capsular_type"], delim=","
        ),
        "qc_checks": functools.partial(
            transpose_multi_columns_parser, transpose=("parameter", ["result", "value", "comment"])
        ),
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
        "amrfinder": functools.partial(parse_amrfinder, amrfinder_uniq_tag=amrfinder_uniq_tag),
        "abricate": functools.partial(
            row_joiner_on_column_parser,
            join_multiple_rows_on=", ",
            columns=["GENE"],
            recode_into_columns=["virulence_database_hits"],
        ),
    }

    print("Processing results...", file=sys.stderr)
    print(f"Results: {results}", file=sys.stderr)
    header_value_dict: dict[str, str] = {"sample": sample_name}
    for result_type, result_path in results.items():
        header_value_dict = header_value_dict | mapping_functions[result_type](result_path)

    if value := aggregate_serotypes(header_value_dict):
        header_value_dict["serotype"] = value

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
        sample_name=snakemake.params.sample_name,
        amrfinder_uniq_tag=snakemake.params.amrfinder_uniq_tag,
    )
