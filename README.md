# Snakemake workflow: `<name>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for `<description>`

After pressing `Use this template`, the following steps are to be performed to finish the configuration of the new repository:

In github:

- Allow Github Actions to create pull requests (in Settings->Actions/general there is a checkbox `Allow GitHub Actions to create and approve pull requests`, check it.)
- Block main branch (in Settings->Branches press `Add branch protection rule`, set `Branch name pattern` to `main` and check the `Require a pull request before merging` checkbox)

In README.md:

- replace `<XY>` values with the correct ones (in badges as well)

## Development

Repository is now created from the template. You can git clone the repository and start the development.

Install snakemake, for example in the environment `snakemake_dev`:

```shell
mamba create -c conda-forge -c bioconda --name snakemake_dev snakemake pre_commit
```

Then set up `pre-commit` in the repository:

```bash
pre-commit install
```

Now before any commit, a defined set of actions will be performed, such as linting, formatting, etc.

### Best practices

Follow the [best practices set by Snakemake authors](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html) ...

### Testing and linting

There is set up an automatic workflow running in the `.tests/` directory. Here you should define a config file and provide any test inputs.
Formatting and linting is also automated for both Python scripts and snakemake rules.

### Conventional PRs

When committing, you must follow the [Conventional Commits spec](https://www.conventionalcommits.org/en/v1.0.0/). Each PR is automatically validated by the GH action.

Further, any push (i.e. after merged PR) to the `main` branch creates in a new PR:

- a new release following the [Semantic Versioning](https://semver.org/)
- an automatic changelog as parsed from the commit history

## Running

Snakemake oworkflows are advised to be run from the workflow directory. Recommended arguments to use, is `--use-conda` to use conda, and the `--conda-prefix` which is the directory where snakemake will install workflow conda environments. Then for example, the command is:

```shell
snakemake --cores {THREADS} --use-conda --conda-prefix {CONDA_DIR}
```

After running, create a [report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html), either `.zip` or `.html`:

```shell
snakemake --report my_first_report.zip
```

It is also advised to [archive](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#sustainable-and-reproducible-archiving) the workflow if successful:

```shell
snakemake --archive my_first_workflow.tar.gz
```

### Recommendations

Use dry-run: `-n` first to confirm that it works as intended:

```shell
snakemake -n
```

To re-run after parameter changes:

```shell
snakemake -n -R `snakemake --list-params-changes`
```

Or after code changes:

```shell
snakemake -n -R `snakemake --list-code-changes`
```

### Debugging

- run with `--notemp` to ignore `temp()` in snakemake rules.
- run with `--show-failed-logs` to automatically show failed log messages

### Issues

- use issues in github to report any problems. This way solving any issue will be stored here for other users or future ourselves.
- Watch the repo to subscribe to pull requests, issues and other repo-related stuff. This way you will be notified of any changes.
