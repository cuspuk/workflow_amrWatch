import os
import sys


def parse_assembly_metrics(bandage_output_file: str):
    dead_ends = None
    contigs = None
    with open(bandage_output_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Dead ends:"):
                dead_ends = int(line.strip().split()[-1])
            elif line.startswith("Node count:"):
                contigs = int(line.strip().split()[-1])
    if not dead_ends or not contigs:
        raise ValueError(f"Could not parse assembly quality from {bandage_output_file=}")
    return dead_ends, contigs


def get_assembly_quality_decision(bandage_output_file: str, max_dead_ends: int, max_contigs: int):
    dead_ends, contigs = parse_assembly_metrics(bandage_output_file)

    if contigs > max_contigs:
        return f"FAIL: Number of contigs is {contigs} which is greater than threshold {max_contigs}"
    if dead_ends > max_dead_ends:
        return f"WARN: Number of dead ends is {dead_ends} which is greater than threshold {max_dead_ends}"
    return f"PASS: Assembly quality fulfills criteria, number of contigs ({contigs}<={max_contigs}) and dead ends ({dead_ends}<={max_dead_ends})"


def evaluate_assembly_quality(bandage_output_file: str, output_path: str, max_dead_ends: int, max_contigs: int):
    decision = get_assembly_quality_decision(bandage_output_file, max_dead_ends, max_contigs)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=False)

    with open(output_path, "w") as f:
        f.write(decision)
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_assembly_quality(
        snakemake.input[0],
        snakemake.output[0],
        snakemake.params.max_dead_ends,
        snakemake.params.max_contigs,
    )
