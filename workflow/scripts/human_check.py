import os
import sys
from dataclasses import dataclass
from enum import StrEnum


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


def get_human_fraction(bracken_file: str, taxonomy_id: str) -> float:
    with open(bracken_file, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            taxa_id, fraction = [line.strip().split("\t")[i] for i in (1, -1)]
            if taxa_id == taxonomy_id:
                return float(fraction)
    return 0.0


def get_human_check_decision(bracken_file: str, taxonomy_id: str, threshold_fraction: float) -> QCRow:
    human_fraction = get_human_fraction(bracken_file, taxonomy_id)
    # genera_log = ",".join([f"{g['genus']}:{g['fraction']}" for g in genera_list])
    human_fraction_str = "{:.3f}".format(human_fraction)
    if human_fraction > threshold_fraction:
        return QCRow(
            QCResult.FAIL,
            "human_contamination",
            human_fraction_str,
            f"Found contamination by {taxonomy_id=} as its fraction is higher than given maximum ({human_fraction} > {threshold_fraction}",
        )
    return QCRow(
        QCResult.PASS,
        "human_contamination",
        human_fraction_str,
        f"No contamination by {taxonomy_id=} as its fraction is lower or equal than maximum ( {human_fraction} <= {threshold_fraction})",
    )


def evaluate_human_contamination(bracken_file: str, output_path: str, taxonomy_id: str, human_fraction: float):
    decision = get_human_check_decision(bracken_file, taxonomy_id, human_fraction)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=False)

    with open(output_path, "w") as f:
        f.write(str(decision))
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_human_contamination(
        snakemake.input[0],
        snakemake.output[0],
        snakemake.params.taxonomy_id,
        snakemake.params.human_fraction,
    )
