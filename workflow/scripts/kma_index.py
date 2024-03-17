import os
import sys

from snakemake.shell import shell


def get_suffix() -> str:
    return snakemake.params.suffix


def get_resfinder_paths() -> list[str]:
    return snakemake.output.resfinder_out


def get_pointfinder_paths() -> list[str]:
    return snakemake.output.pointfinder_out


suffix = get_suffix()
resfinder_names = [os.path.basename(value).removesuffix(suffix) for value in get_resfinder_paths()]
pointfinder_names = [os.path.basename(value).removesuffix(suffix) for value in get_pointfinder_paths()]

sys.stderr = open(snakemake.log[0], "w")

for name in resfinder_names:
    print(f"Processing resfinder: {name}", file=sys.stderr)
    shell(
        "kma_index -i {snakemake.input.resfinder_db}/{name}.fsa -o {snakemake.input.resfinder_db}/{name} >> {snakemake.log} 2>&1"
    )

for name in pointfinder_names:
    print(f"Processing pointfinder: {name}", file=sys.stderr)
    shell(
        "kma_index -i {snakemake.input.pointfinder_db}/{name}/{name}.fsa -o {snakemake.input.pointfinder_db}/{name}/{name} >> {snakemake.log} 2>&1"
    )
