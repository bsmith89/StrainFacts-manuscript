#!/usr/bin/env python3

import sys
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import sfacts as sf
from lib.util import info
from lib.pandas_util import idxwhere


if __name__ == "__main__":
    in_path = sys.argv[1]
    cvrg_path = sys.argv[2]
    species_id = int(sys.argv[3])
    threshold_horizontal_coverage = float(sys.argv[4])
    threshold_dissimilarity = float(sys.argv[5])
    sample_map_path = sys.argv[6]
    sample_to_scg_out_path = sys.argv[7]
    scg_out_path = sys.argv[8]

    raw = sf.data.Metagenotypes.load(in_path)
    cvrg = pd.read_table(
        cvrg_path, index_col=["lib_id", "species_id"], squeeze=True
    ).xs(species_id, level="species_id")
    info(f"Unique barcodes before coverage filtering: {raw.sizes['sample']}")

    filt = raw.mlift("sel", sample=idxwhere(cvrg > threshold_horizontal_coverage))
    ngenotypes = filt.sizes["sample"]

    sample_map = (
        pd.read_table(sample_map_path, index_col="barcode").squeeze().loc[filt.sample]
    )
    info(
        "Unique barcodes before dereplication:",
        sample_map.value_counts().sort_values(ascending=False),
        sep="\n",
    )

    if filt.sizes["sample"] < 2:
        agg = pd.Series([0], index=filt.sample, name="clust")
    else:
        agg = pd.Series(
            (
                AgglomerativeClustering(
                    affinity="cosine",
                    distance_threshold=threshold_dissimilarity,
                    n_clusters=None,
                    linkage="complete",
                ).fit_predict(filt.total_counts())
            ),
            filt.sample,
            name="clust",
        )

    meta = agg.to_frame().join(sample_map)
    info(
        "Barcodes after dereplication:",
        meta.groupby("sample_id")
        .apply(lambda x: len(x.clust.unique()))
        .sort_values(ascending=False),
        sep="\n",
    )
    derep = (
        filt.to_series()
        .reset_index()
        .rename(columns={0: "tally"})
        .join(meta, on="sample")
        .groupby(["sample_id", "clust", "position", "allele"])
        .sum()
        .astype(int)
        .reset_index()
        .assign(label=lambda x: x.sample_id + "_" + x.clust.astype(str))
    )
    sample_to_scg = (
        derep[["sample_id", "label"]].drop_duplicates().set_index("label").sample_id
    )
    sample_to_scg.to_csv(sample_to_scg_out_path, sep="\t", header=False)

    scg = sf.data.Metagenotypes(
        derep.set_index(["label", "position", "allele"])
        .rename_axis(index={"label": "sample"})
        .tally.to_xarray()
    )
    scg.dump(scg_out_path)
