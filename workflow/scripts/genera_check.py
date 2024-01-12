import os
import sys
from typing import TypedDict


class GenusFraction(TypedDict):
    genus: str
    fraction: float


def get_abundant_genera_list(bracken_file: str, threshold_fraction: float) -> list[GenusFraction]:
    genera_list: list[GenusFraction] = []
    with open(bracken_file, "r") as f:
        for line in f.readlines():
            name, fraction = [line.strip().split("\t")[i] for i in (0, -1)]
            if float(fraction) > threshold_fraction:
                genera_list.append(GenusFraction(genus=name, fraction=float(fraction)))

    return genera_list


def get_genera_check_decision(bracken_file: str, threshold_fraction: float) -> str:
    genera_list = get_abundant_genera_list(bracken_file, threshold_fraction)
    genera_log = ",".join([f"{g['genus']}:{g['fraction']}" for g in genera_list])
    if genera_count := len(genera_list) > 1:
        return f"FAIL: There have been {genera_count} genera with fraction > {threshold_fraction}: {genera_log}"
    return f"PASS: Only one genus with fraction > {threshold_fraction}: {genera_log}"


def evaluate_foreign_contamination(bracken_file: str, output_path: str, threshold_fraction: float):
    decision = get_genera_check_decision(bracken_file, threshold_fraction)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=False)

    with open(output_path, "w") as f:
        f.write(decision)
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_foreign_contamination(
        snakemake.input[0],
        snakemake.output[0],
        snakemake.params.fraction_threshold,
    )
