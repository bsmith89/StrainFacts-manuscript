#!/usr/bin/env python3

import sys
import pandas as pd
from tqdm import tqdm
import sqlite3
from lib.pandas_util import idxwhere


if __name__ == "__main__":
    db_path = sys.argv[1]
    agg_path = sys.argv[2]
    ref_path = sys.argv[3]
    threshold_horizontal_coverage = float(sys.argv[4])
    threshold_reference_positions = float(sys.argv[5])
    out_path = sys.argv[6]

    con = sqlite3.connect(db_path)
    agg = pd.read_table(agg_path, index_col="species_id", squeeze=True)
    ref = pd.read_table(ref_path, index_col="species_id", squeeze=True)

    out = {}

    for focal_species in tqdm(idxwhere((agg / ref) > threshold_horizontal_coverage)):
        drplt = pd.read_sql(
            f"""
            SELECT
                lib_id
              , species_position
            FROM snp_x_lib
            JOIN lib USING (lib_id)
            WHERE species_id = '{focal_species}'
            AND lib_type = 'droplet'
            AND (reference_tally + alternative_tally) > 0
        """,
            con=con,
        )
        ref_positions = ref[focal_species]
        if ref_positions < threshold_reference_positions:
            continue
        drplt_positions = drplt.groupby("lib_id").apply(len)
        out[focal_species] = drplt_positions / ref_positions

    out = (
        pd.DataFrame(out)
        .stack(dropna=True)
        .to_frame()
        .reset_index()
        .rename(columns={"level_1": "species_id", 0: "horizontal_coverage"})
        .sort_values(["species_id", "horizontal_coverage"], ascending=(True, False))
    )
    out.to_csv(out_path, sep="\t", index=False)
