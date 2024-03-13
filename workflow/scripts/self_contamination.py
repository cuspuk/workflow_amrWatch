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
        return "result\tparameter\tvalue\tcomment"


def get_ambiguous_rows(vcf_file: str) -> int:
    with open(vcf_file, "r") as f:
        lines = f.readlines()
        unmatched_lines = [line for line in lines if not line.startswith("#")]
    return len(unmatched_lines)


def get_decision(vcf_file: str, max_rows: float, check_level: str) -> QCRow:
    ambiguous_rows = get_ambiguous_rows(vcf_file)
    if ambiguous_rows > max_rows:
        res = QCResult.FAIL if check_level == "FAIL" else QCResult.WARN
        return QCRow(
            res,
            "ambiguous_bases",
            ambiguous_rows,
            f"Number of ambiguous rows is higher than the defined threshold ({ambiguous_rows}>{max_rows})",
        )
    return QCRow(
        QCResult.PASS,
        "ambiguous_bases",
        ambiguous_rows,
        f"Number of ambiguous rows fulfills the criteria ({ambiguous_rows}<={max_rows})",
    )


def evaluate_self_contamination(vcf_file: str, output_path: str, max_rows: float, check_level: str):
    decision = get_decision(vcf_file, max_rows, check_level)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=False)

    with open(output_path, "w") as f:
        f.write(str(decision))
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_self_contamination(
        snakemake.input[0], snakemake.output[0], snakemake.params.max_rows, snakemake.params.check_level
    )
