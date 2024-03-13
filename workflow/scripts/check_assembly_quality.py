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


@dataclass
class AssemblyQuality:
    dead_ends: int
    basepairs: int
    contigs: int

    def _evaluate_base_pairs(self, min_length_in_bp: int, max_length_in_bp: int) -> QCRow:
        parameter = "assembly length"
        if self.basepairs < min_length_in_bp:
            return QCRow(
                QCResult.FAIL,
                parameter,
                self.basepairs,
                f"Number of basepairs in assembly is lower than given threshold ({self.basepairs}<{min_length_in_bp})",
            )
        if self.basepairs > max_length_in_bp:
            return QCRow(
                QCResult.FAIL,
                parameter,
                self.basepairs,
                f"Number of basepairs in assembly is higher than given threshold ({self.basepairs}>{max_length_in_bp})",
            )
        return QCRow(
            QCResult.PASS,
            parameter,
            self.basepairs,
            f"Number of basepairs in assembly fulfills the criteria ({min_length_in_bp}<={self.basepairs}<={max_length_in_bp})",
        )

    def _evaluate_contigs(self, max_contigs: int) -> QCRow:
        if self.contigs > max_contigs:
            return QCRow(
                QCResult.FAIL,
                "number_of_contigs",
                self.contigs,
                f"Number of contigs is greater than threshold ({self.contigs}>{max_contigs})",
            )
        return QCRow(
            QCResult.PASS,
            "number_of_contigs",
            self.contigs,
            f"Number of contigs fulfills criteria ({self.contigs}<={max_contigs})",
        )

    def _evaluate_dead_ends(self, max_dead_ends: int):
        if self.dead_ends > max_dead_ends:
            return QCRow(
                QCResult.WARN,
                "number_of_dead_ends",
                self.dead_ends,
                f"Number of dead ends is greater than threshold ({self.dead_ends}>{max_dead_ends})",
            )
        return QCRow(
            QCResult.PASS,
            "number_of_dead_ends",
            self.dead_ends,
            f"Number of dead ends fulfills criteria ({self.dead_ends}<={max_dead_ends})",
        )

    def evaluate(
        self,
        max_dead_ends: int,
        max_contigs: int,
        min_length_in_bp: int,
        max_length_in_bp: int,
    ):
        rows = [
            self._evaluate_base_pairs(min_length_in_bp, max_length_in_bp),
            self._evaluate_contigs(max_contigs),
            self._evaluate_dead_ends(max_dead_ends),
        ]
        return rows


def parse_assembly_metrics(bandage_output_file: str):
    dead_ends = None
    basepairs = None
    with open(bandage_output_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Dead ends:"):
                dead_ends = int(line.strip().split()[-1])
            elif line.startswith("Total length (bp):"):
                basepairs = int(line.strip().split()[-1])
    if dead_ends is None or basepairs is None:
        raise ValueError(f"Could not parse assembly quality from {bandage_output_file=}")
    return dead_ends, basepairs


def parse_seqkit_stats(seqkit_stats_file: str):
    with open(seqkit_stats_file, "r") as f:
        lines = f.readlines()
        header = lines[0].strip().split("\t")
        try:
            num_seqs_index = header.index("num_seqs")
        except ValueError:
            raise ValueError(f"Could not find num_seqs in {seqkit_stats_file=}")
        row = lines[1].strip().split("\t")
        return int(row[num_seqs_index])


def get_assembly_quality_decision(
    bandage_output_file: str,
    seqkit_stats_file: str,
):
    dead_ends, basepairs = parse_assembly_metrics(bandage_output_file)
    contigs = parse_seqkit_stats(seqkit_stats_file)
    return AssemblyQuality(dead_ends, basepairs, contigs)


def evaluate_assembly_quality(
    bandage_output_file: str,
    seqkit_stats_file: str,
    output_path: str,
    max_dead_ends: int,
    max_contigs: int,
    min_length_in_bp: int,
    max_length_in_bp: int,
):
    parsed_values = get_assembly_quality_decision(bandage_output_file, seqkit_stats_file)
    decisions = parsed_values.evaluate(max_dead_ends, max_contigs, min_length_in_bp, max_length_in_bp)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=True)

    with open(output_path, "w") as f:
        for decision in decisions:
            f.write(str(decision))
            f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_assembly_quality(
        snakemake.input.tsv,
        snakemake.input.stats,
        snakemake.output[0],
        snakemake.params.max_dead_ends,
        snakemake.params.max_contigs,
        snakemake.params.min_length_in_bp,
        snakemake.params.max_length_in_bp,
    )
