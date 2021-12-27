#!/usr/bin/env python3

import sys
import pandas as pd
import sfacts as sf
from lib.pandas_util import idxwhere

if __name__ == "__main__":
    mgen_path = sys.argv[1]
    scg_path = sys.argv[2]
    fit_path = sys.argv[3]
    scg_to_sample_path = sys.argv[4]
    library_to_sample_path = sys.argv[5]
    rabund_threshold = float(sys.argv[6])
    out_path = sys.argv[7]

    mgen = sf.data.Metagenotypes.load(mgen_path)
    drplt = sf.data.Metagenotypes.load(scg_path)
    inference = sf.data.World.load(fit_path)

    mgen_consensus = mgen.to_estimated_genotypes()
    scg = drplt.to_estimated_genotypes()
    inferred_geno = inference.genotypes
    inferred_comm = inference.communities

    position = list(set(scg.position.values) & set(inferred_geno.position.values))
    mgen_consensus, scg, inferred_geno = [
        x.mlift("sel", position=position) for x in [mgen_consensus, scg, inferred_geno]
    ]

    mgen_entropy = mgen.entropy("sample")
    scg_entropy = drplt.entropy("sample")
    comm_entropy = inferred_comm.entropy("sample")
    geno_entropy = inferred_geno.entropy("strain")
    scg_horizontal_coverage = drplt.horizontal_coverage()
    mgen_horizontal_coverage = mgen.horizontal_coverage()

    library_to_sample = pd.read_table(
        library_to_sample_path, index_col="barcode", squeeze=True
    )
    scg_to_sample = pd.read_table(
        scg_to_sample_path, names=["strain", "sample"], index_col="strain", squeeze=True
    )
    mgen_to_sample = library_to_sample.loc[mgen.sample]


    fdist_any_strain = sf.match_genotypes(
        scg.to_world(), inferred_geno.to_world()
    )[1]
    ddist_any_strain = sf.match_genotypes(
        scg.to_world(), inferred_geno.discretized().to_world()
    )[1]
    fdist_any_mgen = sf.match_genotypes(
        scg.to_world(), mgen_consensus.to_world()
    )[1]
    ddist_any_mgen = sf.match_genotypes(
        scg.to_world(), mgen_consensus.discretized().to_world()
    )[1]

    out = []
    for focal_mgen in idxwhere(mgen_to_sample.isin(scg_to_sample.unique())):
        focal_sample = mgen_to_sample[focal_mgen]
        focal_strains = idxwhere((inferred_comm.sel(sample=focal_mgen) > rabund_threshold).to_series())
        fdist_focal_strain = sf.match_genotypes(
            scg.to_world(), inferred_geno.mlift('sel', strain=focal_strains).to_world()
        )[1]
        ddist_focal_strain = sf.match_genotypes(
            scg.to_world(), inferred_geno.mlift('sel', strain=focal_strains).discretized().to_world()
        )[1]
        fdist_focal_mgen = sf.match_genotypes(
            scg.to_world(), mgen_consensus.mlift('sel', strain=[focal_mgen]).to_world()
        )[1]
        ddist_focal_mgen = sf.match_genotypes(
            scg.to_world(), mgen_consensus.mlift('sel', strain=[focal_mgen]).discretized().to_world()
        )[1]
        for focal_scg in idxwhere(scg_to_sample == focal_sample):
            out.append(dict(
                scg=focal_scg,
                sample=focal_sample,
                mgen=focal_mgen,
                mgen_entropy=mgen_entropy.sel(sample=focal_mgen).values,
                mgen_horizontal_coverage=mgen_horizontal_coverage.sel(sample=focal_mgen).values,
                scg_entropy=scg_entropy.sel(sample=focal_scg).values,
                scg_horizontal_coverage=scg_horizontal_coverage.sel(sample=focal_scg).values,
                comm_entropy=comm_entropy.sel(sample=focal_mgen).values,
                fdist_any_strain=fdist_any_strain[focal_scg],
                ddist_any_strain=ddist_any_strain[focal_scg],
                fdist_any_mgen=fdist_any_mgen[focal_scg],
                ddist_any_mgen=ddist_any_mgen[focal_scg],
                fdist_focal_strain=fdist_focal_strain[focal_scg],
                ddist_focal_strain=ddist_focal_strain[focal_scg],
                fdist_focal_mgen=fdist_focal_mgen[focal_scg],
                ddist_focal_mgen=ddist_focal_mgen[focal_scg],
        ))
    out = pd.DataFrame(out)

    out.to_csv(out_path, sep='\t', index=False)
