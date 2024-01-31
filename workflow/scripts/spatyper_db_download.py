import os
import sys

import spaTyper


def download_db_for_spaTyper(db_dir: str):
    os.makedirs(db_dir, exist_ok=True)

    spaTyper.utils.download_file_repeats(db_dir, False)
    spaTyper.utils.download_file_types(db_dir, False)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr
    download_db_for_spaTyper(
        snakemake.params.db_dir,
    )
