import os
import sys


def parse_assembly_metrics(bandage_output_file: str):
    dead_ends = None
    contigs = None
    basepairs = None
    with open(bandage_output_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Dead ends:"):
                dead_ends = int(line.strip().split()[-1])
            elif line.startswith("Node count:"):
                contigs = int(line.strip().split()[-1])
            elif line.startswith("Total length (bp):"):
                basepairs = int(line.strip().split()[-1])
    if dead_ends is None or contigs is None or basepairs is None:
        raise ValueError(f"Could not parse assembly quality from {bandage_output_file=}")
    return dead_ends, contigs, basepairs


def get_assembly_quality_decision(
    bandage_output_file: str, max_dead_ends: int, max_contigs: int, min_length_in_bp: int, max_length_in_bp: int
):
    dead_ends, contigs, basepairs = parse_assembly_metrics(bandage_output_file)

    if basepairs < min_length_in_bp:
        return f"FAIL: Number of basepairs in assembly is lower than given threshold ({basepairs}<{min_length_in_bp})"
    if basepairs > max_length_in_bp:
        return f"FAIL: Number of basepairs in assembly is higher than given threshold ({basepairs}>{max_length_in_bp})"
    if contigs > max_contigs:
        return f"FAIL: Number of contigs is {contigs} which is greater than threshold {max_contigs}"
    if dead_ends > max_dead_ends:
        return f"WARN: Number of dead ends is {dead_ends} which is greater than threshold {max_dead_ends}"
    return f"PASS: Assembly quality fulfills criteria, assembly length ({max_length_in_bp}>={basepairs}>={min_length_in_bp}), number of contigs ({contigs}<={max_contigs}) and dead ends ({dead_ends}<={max_dead_ends})"


def evaluate_assembly_quality(
    bandage_output_file: str,
    output_path: str,
    max_dead_ends: int,
    max_contigs: int,
    min_length_in_bp: int,
    max_length_in_bp: int,
):
    decision = get_assembly_quality_decision(
        bandage_output_file, max_dead_ends, max_contigs, min_length_in_bp, max_length_in_bp
    )

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=True)

    with open(output_path, "w") as f:
        f.write(decision)
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_assembly_quality(
        snakemake.input.tsv,
        snakemake.output[0],
        snakemake.params.max_dead_ends,
        snakemake.params.max_contigs,
        snakemake.params.min_length_in_bp,
        snakemake.params.max_length_in_bp,
    )
