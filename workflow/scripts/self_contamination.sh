#!/usr/bin/env bash


rows=$(grep -cv '^#' "${snakemake_input[0]}" 2> "${snakemake_log[0]}" );


out_dir="$(dirname "${snakemake_output[0]}" 2>> "${snakemake_log[0]}")"

mkdir -p "$out_dir" 2>> "${snakemake_log[0]}";

if [ "$rows" -gt "${snakemake_params[max_rows]}" ]
then
    echo "PASS: Number of ambiguous rows fulfills the criteria ($rows<=${snakemake_params[max_rows]})" > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"
else
    echo "FAIL: Number of ambiguous rows is higher than the defined threshold ($rows>${snakemake_params[max_rows]})" > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"
fi
