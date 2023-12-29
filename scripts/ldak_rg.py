import concurrent.futures
import pathlib
import shlex
import subprocess


def main():
    root = pathlib.Path("/data1/home/mnz2108/git/maxgcp-analysis/data/gwas_results")
    files = list(root.glob("plink_white_british*.glm.linear.summaries"))

    rows = [
        f"ldak_rg {file_1.as_posix()} {file_2.as_posix()}\n"
        for file_1, file_2 in itertools.combinations(files, 2)
    ]
    with open("rg_jobs.txt", "w") as f:
        f.writelines(rows)
