import os
import shutil
from tempfile import TemporaryDirectory

from snakemake.shell import shell

if os.path.exists(snakemake.params.dir):
    shutil.rmtree(snakemake.params.dir)

os.makedirs(snakemake.params.dir, exist_ok=True)
os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory() as tempdir:
    tmp_zip = os.path.join(tempdir, "test.zip")
    tmp_unzipped = os.path.join(tempdir, "tmp_unzipped")
    shell(
        "(mkdir -p {tmp_unzipped} &&"
        " curl -SL {snakemake.params.repo} -o {tmp_zip} &&"
        " unzip {tmp_zip} -d {tmp_unzipped} &&"
        " mv {tmp_unzipped}/{snakemake.params.name_after_unzip}/* {snakemake.params.dir}"
        " ) {log}"
    )
