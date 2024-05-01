import pathlib

import polars as pl


def main():
    root = pathlib.Path("~/git/maxgcp-analysis/data/phenotypes")

    (
        pl.read_csv(root / "pheno.tsv", separator="\t")
        .select("#FID", "IID", pl.col("^b_.+$").add(2))
        .write_csv(root / "binary_pheno_recoded.tsv", separator="\t")
    )


if __name__ == "__main__":
    main()
