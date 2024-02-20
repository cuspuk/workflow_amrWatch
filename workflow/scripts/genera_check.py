import os
import sys
from dataclasses import dataclass
from enum import StrEnum
from typing import TypedDict


class QCResult(StrEnum):
    PASS = "PASS"
    FAIL = "FAIL"
    WARN = "WARN"


@dataclass
class QCRow:
    result: QCResult
    parameter: str
    value: int | float | str
    comment: str

    def __str__(self):
        return f"{self.result}\t{self.parameter}\t{self.value}\t{self.comment}"

    @staticmethod
    def header():
        return f"result\tparameter\tvalue\tcomment"


class GenusFraction(TypedDict):
    genus: str
    fraction: float


def get_abundant_genera_list(bracken_file: str, threshold_fraction: float) -> list[GenusFraction]:
    genera_list: list[GenusFraction] = []
    with open(bracken_file, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            name, fraction = [line.strip().split("\t")[i] for i in (0, -1)]
            if float(fraction) > threshold_fraction:
                genera_list.append(GenusFraction(genus=name, fraction=float(fraction)))

    return genera_list


def get_genera_check_decision(bracken_file: str, threshold_fraction: float) -> QCRow:
    genera_list = get_abundant_genera_list(bracken_file, threshold_fraction)
    genera_log = ",".join([f"{g['genus']}:{g['fraction']}" for g in genera_list])
    if (genera_count := len(genera_list)) > 1:
        return QCRow(
            QCResult.FAIL,
            "foreign_contamination",
            genera_log,
            f"There have been {genera_count} genera with fraction > {threshold_fraction}",
        )
    return QCRow(
        QCResult.PASS, "foreign_contamination", genera_log, f"Only one genus with fraction > {threshold_fraction}"
    )


def evaluate_foreign_contamination(bracken_file: str, output_path: str, threshold_fraction: float):
    decision = get_genera_check_decision(bracken_file, threshold_fraction)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=False)

    with open(output_path, "w") as f:
        f.write(str(decision))
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_foreign_contamination(
        snakemake.input[0],
        snakemake.output[0],
        snakemake.params.fraction_threshold,
    )
