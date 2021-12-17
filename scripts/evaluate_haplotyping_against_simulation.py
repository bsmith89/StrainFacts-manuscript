#!/usr/bin/env python3

import sfacts as sf
import sys
import pandas as pd


if __name__ == "__main__":
    sim_path = sys.argv[1]
    fit_path = sys.argv[2]
    bench_path = sys.argv[3]
    out_path = sys.argv[4]
    fit = sf.data.World.load(fit_path)
    sim = sf.data.World.open(sim_path).sel(position=fit.position, sample=fit.sample)
    bench = pd.read_table(bench_path)

    out = dict(
        metagenotype_prediction_error=sf.evaluation.metagenotype_error2(fit)[0],
        metagenotype_prediction_discretized_error=sf.evaluation.metagenotype_error2(fit, discretized=True)[0],
        braycurtis_trans_error=sf.evaluation.braycurtis_error(sim, fit)[0],
        unifrac_trans_error=sf.evaluation.unifrac_error2(sim, fit)[0],
        unifrac_cis_error=sf.evaluation.unifrac_error(sim, fit)[0],
        rank_abundance_error=sf.evaluation.rank_abundance_error(sim, fit)[0],
        community_entropy_error=sf.evaluation.community_entropy_error(sim, fit)[0],
        fwd_genotype_error=sf.evaluation.weighted_genotype_error(sim, fit),
        rev_genotype_error=sf.evaluation.weighted_genotype_error(fit, sim),
        fwd_discrete_genotype_error=sf.evaluation.discretized_weighted_genotype_error(sim, fit),
        rev_discrete_genotype_error=sf.evaluation.discretized_weighted_genotype_error(fit, sim),
        runtime=bench["s"][0],
    )
    (
        pd.DataFrame(out, index=[fit_path])
        .rename_axis(index="fit_path")
        .to_csv(out_path, sep="\t")
    )
