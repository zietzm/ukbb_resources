import polars as pl


def main():
    root = "~/git/maxgcp-analysis/data/phenotypes"
    path = f"{root}/pheno_regenie.tsv"
    df = pl.read_csv(path, separator="\t")

    # Remove phenotypes with no variance
    columns_to_keep = ["FID", "IID"] + [
        col for col in df.select(pl.col("^(b|q)_.+$")).columns if df[col].std() != 0
    ]

    df = df.select(columns_to_keep)
    out_path = f"{root}/pheno_regenie_filtered.tsv"
    df.write_csv(out_path, separator="\t", null_value="NA")


if __name__ == "__main__":
    main()
