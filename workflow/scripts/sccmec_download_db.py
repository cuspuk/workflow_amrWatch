import os
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if not os.path.exists(snakemake.params.db_dir):
    os.makedirs(snakemake.params.db_dir, exist_ok=True)

with TemporaryDirectory() as tempdir:
    tmp_zip = os.path.join(tempdir, "test.zip")
    tmp_unzipped = os.path.join(tempdir, "tmp_unzipped")
    shell(
        "(mkdir -p {tmp_unzipped} &&"
        " curl -SL {snakemake.params.repo} -o {tmp_zip} &&"
        " unzip {tmp_zip} -d {tmp_unzipped} &&"
        " mv {tmp_unzipped}/{snakemake.params.repo_db_name}/* {snakemake.output.db_dir}"
        " ) {log}"
    )
