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
resfinder_out_base = get_resfinder_paths()[0].removesuffix(suffix)
pointfinder_out_base = get_pointfinder_paths()[0].removesuffix(suffix)

os.makedirs(resfinder_out_base, exist_ok=True)
os.makedirs(pointfinder_out_base, exist_ok=True)

sys.stderr = open(snakemake.log[0], "w")

for name in resfinder_names:
    print(f"Processing resfinder: {name}", file=sys.stderr)
    out_name = os.path.join(resfinder_out_base, name)
    shell("kma_index -i {snakemake.input.resfinder_db}/{name}.fsa -o {out_name} >> {snakemake.log} 2>>&1")

for name in pointfinder_names:
    print(f"Processing pointfinder: {name}", file=sys.stderr)
    out_name = os.path.join(pointfinder_out_base, name)
    shell("kma_index -i {snakemake.input.pointfinder_db}/{name}.fsa -o {out_name} >> {snakemake.log} 2>>&1")
