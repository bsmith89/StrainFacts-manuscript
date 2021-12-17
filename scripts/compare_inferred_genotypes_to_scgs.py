#!/usr/bin/env python3

import sys
import pandas as pd
import sfacts as sf
import numpy as np
from lib.pandas_util import idxwhere


if __name__ == "__main__":
    mgen_path = sys.argv[1]
    scg_path = sys.argv[2]
    fit_path = sys.argv[3]
    scg_to_sample_path = sys.argv[4]
    library_to_sample_path = sys.argv[5]
    threshold = float(sys.argv[6])
    pseudo = float(sys.argv[7])  # e.g. 1e-10
    out_path = sys.argv[8]

    mgen = sf.data.Metagenotypes.load(mgen_path)
    drplt = sf.data.Metagenotypes.load(scg_path)
    inference = sf.data.World.load(fit_path)

    mgen_consensus = mgen.to_estimated_genotypes(pseudo=pseudo)
    scg = drplt.to_estimated_genotypes(pseudo=pseudo)
    inferred_geno = inference.genotypes
    inferred_comm = inference.communities

    position = list(set(scg.position.values) & set(inferred_geno.position.values))
    mgen_consensus, scg, inferred_geno = [
        x.mlift("sel", position=position) for x in [mgen_consensus, scg, inferred_geno]
    ]

    # scg_to_sample = scg.strain.to_series().str.split("_").apply(lambda x: x[0]) + ".m"
    library_to_sample = pd.read_table(
        library_to_sample_path, index_col="barcode", squeeze=True
    )
    scg_to_sample = pd.read_table(
        scg_to_sample_path, names=["strain", "sample"], index_col="strain", squeeze=True
    )

    smallest_fdist_all_strains = sf.match_genotypes(
        scg.to_world(), inferred_geno.to_world()
    )[1]
    smallest_ddist_all_strains = sf.match_genotypes(
        scg.to_world(), inferred_geno.discretized().to_world()
    )[1]
    smallest_fdist_all_mgen = sf.match_genotypes(
        scg.to_world(), mgen_consensus.to_world()
    )[1]
    smallest_ddist_all_mgen = sf.match_genotypes(
        scg.to_world(), mgen_consensus.discretized().to_world()
    )[1]

    scg_horizontal_coverage = (drplt.total_counts() > 0).mean("position")

    focal_samples = []
    smallest_fdist_focal_strains = []
    smallest_ddist_focal_strains = []
    smallest_fdist_focal_mgen = []
    smallest_ddist_focal_mgen = []
    horizontal_coverage_focal_mgen = []
    community_entropy_focal_sample = []

    for focal_sample, d in scg_to_sample.reset_index(name="sample").groupby("sample"):
        focal_sample_library_list = idxwhere(
            (library_to_sample == focal_sample)
            & (library_to_sample.index.to_series().isin(inferred_comm.sample.values))
        )
        if not focal_sample_library_list:
            continue
        focal_sample_scg_list = scg_to_sample[scg_to_sample == focal_sample].index
        focal_sample_strain_list = idxwhere(
            (inferred_comm.sel(sample=focal_sample_library_list) > threshold)
            .any("sample")
            .to_series()
        )

        focal_samples.append(pd.Series([focal_sample]*len(d), index=focal_sample_scg_list))
        smallest_fdist_focal_strains.append(
            sf.match_genotypes(
                scg.mlift("sel", strain=focal_sample_scg_list).to_world(),
                inferred_geno.mlift("sel", strain=focal_sample_strain_list).to_world(),
            )[1]
        )
        smallest_ddist_focal_strains.append(
            sf.match_genotypes(
                scg.mlift("sel", strain=focal_sample_scg_list).to_world(),
                inferred_geno.mlift("sel", strain=focal_sample_strain_list)
                .discretized()
                .to_world(),
            )[1]
        )
        smallest_fdist_focal_mgen.append(
            sf.match_genotypes(
                scg.mlift("sel", strain=focal_sample_scg_list).to_world(),
                mgen_consensus.mlift(
                    "sel", strain=focal_sample_library_list
                ).to_world(),
            )[1]
        )
        smallest_ddist_focal_mgen.append(
            sf.match_genotypes(
                scg.mlift("sel", strain=focal_sample_scg_list).to_world(),
                mgen_consensus.mlift("sel", strain=focal_sample_library_list)
                .discretized()
                .to_world(),
            )[1]
        )
        horizontal_coverage_focal_mgen.append(
            pd.Series(
                float(
                    (
                        mgen.total_counts()
                        .sel(sample=focal_sample_library_list)
                        .sum("sample")
                        > 0
                    ).mean()
                ),
                index=d.strain,
            )
        )
        community_entropy_focal_sample.append(
            pd.Series(
                float(
                    (
                        inferred_comm.entropy("sample")
                        .sel(sample=focal_sample_library_list)
                        .mean(
                            "sample"
                        )  # FIXME: The mean entropy across multiple replicate libraries isn't really what I care about ...
                    )
                ),
                index=d.strain,
            )
        )

    # breakpoint()
    if focal_samples:
        focal_sample = pd.concat(focal_samples)
    else:
        focal_sample = np.nan

    if smallest_fdist_focal_strains:
        smallest_fdist_focal_strains = pd.concat(smallest_fdist_focal_strains)
    else:
        smallest_fdist_focal_strains = np.nan

    if smallest_ddist_focal_strains:
        smallest_ddist_focal_strains = pd.concat(smallest_ddist_focal_strains)
    else:
        smallest_ddist_focal_strains = np.nan

    if smallest_fdist_focal_mgen:
        smallest_fdist_focal_mgen = pd.concat(smallest_fdist_focal_mgen)
    else:
        smallest_fdist_focal_mgen = np.nan

    if smallest_ddist_focal_mgen:
        smallest_ddist_focal_mgen = pd.concat(smallest_ddist_focal_mgen)
    else:
        smallest_ddist_focal_mgen = np.nan

    if horizontal_coverage_focal_mgen:
        horizontal_coverage_focal_mgen = pd.concat(horizontal_coverage_focal_mgen)
    else:
        horizontal_coverage_focal_mgen = np.nan

    if community_entropy_focal_sample:
        community_entropy_focal_sample = pd.concat(community_entropy_focal_sample)
    else:
        community_entropy_focal_sample = np.nan

    # breakpoint()
    out = pd.DataFrame(
        dict(
            focal_sample=focal_sample,
            smallest_fdist_all_strains=smallest_fdist_all_strains,
            smallest_ddist_all_strains=smallest_ddist_all_strains,
            smallest_fdist_all_mgen=smallest_fdist_all_mgen,
            smallest_ddist_all_mgen=smallest_ddist_all_mgen,
            smallest_fdist_focal_strains=smallest_fdist_focal_strains,
            smallest_ddist_focal_strains=smallest_ddist_focal_strains,
            smallest_fdist_focal_mgen=smallest_fdist_focal_mgen,
            smallest_ddist_focal_mgen=smallest_ddist_focal_mgen,
            scg_horizontal_coverage=scg_horizontal_coverage,
            horizontal_coverage_focal_mgen=horizontal_coverage_focal_mgen,
            community_entropy_focal_sample=community_entropy_focal_sample,
        )
    )
    out.rename_axis(index="scg").to_csv(out_path, sep="\t")
