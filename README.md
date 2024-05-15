# Snakemake workflow: amrWatch

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/cuspuk/workflow_amrWatch/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/cuspuk/workflow_amrWatch/actions?query=branch%3Amain+workflow%3ATests)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

A Snakemake workflow for antimicrobial resistance samples.

## Installing and running

To install the workflow, simply git clone the repository.

Then, install the following conda environment:

```bash
mamba create -c conda-forge -c bioconda --name snakemake_v7_25 python=3.11 snakemake=7.25 peppy snakemake-wrapper-utils
```

**IMPORTANT**: change the directory to the cloned repository - workflow directory. Every relative path mentioned is relative to this directory.

### Preparing data and configuration

First, prepare your data configuration using PEP file, see [samples.csv](config/pep/samples.csv) as an example. Create new file or update the existing one.

Second, you must configure the `config/config.yaml` file:

- Do not forget to update the `pepfile` path, if you created your own PEP file.
- Provide the path to reference input data (at the config header).
- Update other values as desired.

### Running the workflow

There are many arguments to use when running a snakemake workflow, see [the documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html). Recommended arguments to use, is `--use-conda` to use conda, and the `--conda-prefix` which is the directory where snakemake will install workflow conda environments. Then for example, the command is:

First, it is advised to dry-run the snakemake workflow using `--dry-run`, i.e.:

```shell
snakemake --cores {THREADS} --use-conda --rerun-incomplete --printshellcmds --dry-run
```

A basic summary of outputs and rules is outputted for you to verify. Then, run snakemake without `--dry-run`

```shell
snakemake --cores {THREADS} --use-conda --rerun-incomplete --printshellcmds
```

### Debugging

- run with `--notemp` to ignore `temp()` in snakemake rules.
- run with `--show-failed-logs` to automatically show failed log messages.

### Issues

- use issues in github to report any problems. This way solving any issue will be stored here for other users or future ourselves.
- Watch the repo to subscribe to pull requests, issues and other repo-related stuff. This way you will be notified of any changes.

## Development

Install `snakemake` with `pre_commit`, for example in the environment `snakemake_dev`:

```shell
mamba create -c conda-forge -c bioconda --name snakemake_dev snakemake=7.25 snakemake-wrapper-utils pre_commit peppy
```

Then set up `pre-commit` in the cloned repository:

```bash
pre-commit install
```

Now before any commit, a defined set of actions will be performed, such as linting, formatting, etc.

### Testing and linting

There is set up an automatic workflow running in the `.tests/` directory. Here you should define a config file and provide any test inputs.
Formatting and linting is also automated for both Python scripts and snakemake rules.

### Conventional PRs

When committing, you must follow the [Conventional Commits spec](https://www.conventionalcommits.org/en/v1.0.0/). Each PR is automatically validated by the GH action.

Further, any push (i.e. after merged PR) to the `main` branch creates in a new PR:

- a new release following the [Semantic Versioning](https://semver.org/)
- an automatic changelog as parsed from the commit history
