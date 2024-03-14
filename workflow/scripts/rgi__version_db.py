import json
import os
import sys

# Script to parse RGI card_json and store version to a file
# Example: {"card_json": {"data_version": "3.2.8", "model_type_used": ["homolog", "variant", "rRNA", "overexpression", "knockout"]}, "card_canonical": {"data_version": "N/A", "model_type_used": []}, "card_variants": {"data_version": "N/A", "model_type_used": []}, "card_kmers": {"kmer_sizes": []}}


def parse_rgi_version(card_json_path: str):
    print("Parsing RGI version from", card_json_path, file=sys.stderr)
    with open(card_json_path, "r") as f:
        json_dict = json.load(f)
    print("Parsed RGI version:", json_dict["card_json"]["data_version"], file=sys.stderr)
    return json_dict["card_json"]["data_version"]


def store_version(version: str, output_path: str):
    print("Storing RGI version to", output_path, file=sys.stderr)
    if not os.path.isdir(os.path.dirname(output_path)):
        print("Creating directory", os.path.dirname(output_path), file=sys.stderr)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as f:
        f.write(version)


def get_rgi_version(card_json_path: str, output_path: str):
    version = parse_rgi_version(card_json_path)
    store_version(version, output_path)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    get_rgi_version(
        snakemake.params.json,
        snakemake.output[0],
    )
