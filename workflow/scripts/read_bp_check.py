import os
import sys
from dataclasses import dataclass
from enum import StrEnum

CUTADAPT_LINE_PREFIX = "        Total written (filtered):"
CUTADAPT_UNITS = "bp"


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


def extract_number_of_basepairs(line: str):
    line = line.removeprefix(CUTADAPT_LINE_PREFIX)
    line = line.split(CUTADAPT_UNITS)[0]
    return int(line.replace(",", ""))


def is_line_correct(line: str):
    return line.startswith(CUTADAPT_LINE_PREFIX) and CUTADAPT_UNITS in line


def get_decision_for_basepairs(filename: str, min_bp: int):
    with open(filename, "r") as file_one:
        for line in file_one:
            if is_line_correct(line):
                bp = extract_number_of_basepairs(line)
                if bp < min_bp:
                    return QCRow(
                        QCResult.FAIL, "min_bp", bp, f"Number of basepairs is less than threshold ({bp}<{min_bp})"
                    )
                else:
                    return QCRow(
                        QCResult.PASS, "min_bp", bp, f"Number of basepairs fulfills the criteria ({bp}>={min_bp})"
                    )
    raise ValueError(f"Regex={CUTADAPT_LINE_PREFIX} did not matched any row in cutadapt file={filename}")


def evaluate_minbp_check(filename: str, output_path: str, min_bp: int):
    decision = get_decision_for_basepairs(filename, min_bp)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=False)

    with open(output_path, "w") as f:
        f.write(str(decision))
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_minbp_check(
        snakemake.input[0],
        snakemake.output[0],
        snakemake.params.min_bp,
    )
