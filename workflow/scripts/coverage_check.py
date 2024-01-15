import locale
import os
import sys

_REGEX_VALUE = "     mean coverageData = "


def convert_qualimap_coverage_value(string_value: str):
    locale.setlocale(locale.LC_ALL, "en_US.UTF-8")
    return locale.atof(string_value.removesuffix("X"))


def found_qualimap_row(line: str) -> bool:
    if line.startswith(_REGEX_VALUE):
        return True
    return False


def clean_qualimap_row(line: str) -> str:
    return line.removeprefix(_REGEX_VALUE).split()[0]


def parse_mean_coverage(line: str) -> float:
    float_str = clean_qualimap_row(line)
    return convert_qualimap_coverage_value(float_str)


def parse_qualimap_for_coverage(filename: str) -> float:
    with open(filename, "r") as file_one:
        for line in file_one:
            if found_qualimap_row(line):
                return parse_mean_coverage(line)
    raise ValueError(f"Regex={_REGEX_VALUE} did not matched any row in qualimap file={filename}")


def get_decision_for_coverage(coverage: float, warn_threshold: float, fail_threshold: float) -> str:
    if coverage < fail_threshold:
        return f"FAIL: Mean {coverage=} is lower than minimum threshold={fail_threshold}"
    if coverage < warn_threshold:
        return f"WARN: Mean {coverage=} is higher than minimum threshold={fail_threshold} but lower than good threshold={warn_threshold}"

    return f"PASS: Mean {coverage=} is higher than good threshold={warn_threshold}"


def evaluate_coverage_check(filename: str, output_path: str, warn_threshold: float, fail_threshold: float):
    decision = get_decision_for_coverage(parse_qualimap_for_coverage(filename), warn_threshold, fail_threshold)

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o777, exist_ok=False)

    with open(output_path, "w") as f:
        f.write(decision)
        f.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_coverage_check(
        snakemake.params.genome_results_file,
        snakemake.output[0],
        snakemake.params.warn_threshold,
        snakemake.params.fail_threshold,
    )
