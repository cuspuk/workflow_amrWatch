import os
import sys


def get_mlst(mlst_file: str, mlst_index: int) -> str:
    with open(mlst_file, "r") as f:
        row = f.readline().rstrip().split("\t")
        return row[mlst_index]


def get_clonal_complex(mlst_value: str, profile_tsv_file: str) -> tuple[list[str], list[str]]:
    with open(profile_tsv_file, "r") as f:
        header = f.readline().rstrip().split("\t")
        st_index = header.index("ST")
        rows = [row.rstrip().split("\t") for row in f.readlines()]
    for row in rows:
        if row[st_index] == mlst_value:
            return header, row
    return (["ST", "clonal_complex"], [mlst_value, "ST_not_found"])


def get_and_output(mlst_file: str, mlst_index: int, profile_tsv_file: str, output: str):
    mlst_value = get_mlst(mlst_file, mlst_index)
    header, row = get_clonal_complex(mlst_value, profile_tsv_file)

    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    get_and_output(
        snakemake.input.mlst,
        snakemake.params.mlst_index,
        snakemake.input.profile_tsv,
        snakemake.output[0],
    )
