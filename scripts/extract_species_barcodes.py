#!/usr/bin/env python3
import pandas as pd
import sqlite3
import sys
import sfacts as sf


def metagenotype_db_to_xarray(df):
    """Convert from sc-validate-haplotypes project database schema to a metagenotype."""
    return (
        df.rename_axis(columns="allele")
        .rename(columns=dict(alternative_tally="alt", reference_tally="ref"))
        .rename_axis(index=dict(lib_id="sample", species_position="position"))
        .stack()
        .to_xarray()
        .fillna(0)
        .sortby("allele")
    )


if __name__ == "__main__":
    db_path = sys.argv[1]
    species_id = sys.argv[2]
    lib_type = sys.argv[3]
    out_path = sys.argv[4]

    con = sqlite3.connect(db_path)
    data = pd.read_sql(
        f"""
SELECT
  lib_id,
  species_position,
  reference_tally,
  alternative_tally
FROM snp_x_lib
JOIN lib USING (lib_id)
WHERE species_id = '{species_id}'
  AND lib_type = '{lib_type}'
    """,
        con=con,
        index_col=["lib_id", "species_position"],
    )
    mgen = sf.data.Metagenotypes(metagenotype_db_to_xarray(data))
    mgen.dump(out_path)
