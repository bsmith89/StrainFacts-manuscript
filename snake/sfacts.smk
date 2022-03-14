import math


use rule start_jupyter as start_jupyter_sfacts with:
    conda:
        "conda/sfacts.yaml"


use rule start_ipython as start_ipython_sfacts with:
    conda:
        "conda/sfacts.yaml"


use rule start_shell as start_shell_sfacts with:
    conda:
        "conda/sfacts.yaml"


rule compile_species_variation_from_vcf:
    output:
        "data/gtprodb.sp-{species}.genotype.nc",
    input:
        script="scripts/vcf_to_gtpro_pileup.py",
        gtpro_snp_dict="ref/gtpro/variants_main.covered.hq.snp_dict.tsv",
        vcf="raw/gtpro_refs/variation_in_species/{species}/core_snps.vcf.gz",
    shell:
        """
        {input.script} {input.gtpro_snp_dict} {input.vcf} {wildcards.species} {output}
        """


rule extract_species_barcodes:
    output:
        "data/ucfmt.sp-{species_id}.{lib_type}.nc",
    input:
        script="scripts/extract_species_barcodes.py",
        db="data/ucfmt.db",
    wildcard_constraints:
        lib_type="metagenotype|genotype",
        species_id=noperiod_wc,
    params:
        lib_type=lambda w: {"genotype": "droplet", "metagenotype": "metagenome"}[
            w.lib_type
        ],
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.db} {wildcards.species_id} {params.lib_type} {output}
        """


rule filter_and_dereplicate_barcodes_by_coverage:
    output:
        scg="{stem}.sp-{species_id}.derep.genotype.nc",
        sample_to_scg="{stem}.sp-{species_id}.derep.barcode_to_sample.tsv",
    input:
        script="scripts/dereplicate_barcodes_by_coverage.py",
        drplt="{stem}.sp-{species_id}.genotype.nc",
        cvrg="{stem}.genotype.horizontal_coverage.tsv",
        sample_map="data/ucfmt.barcode_to_sample.tsv",
    params:
        species_id=lambda w: w.species_id,
        threshold_horizontal_coverage=0.01,
        threshold_dissimilarity=0.5,
    conda:
        "conda/sfacts.yaml"
    shell:
        r"""
        {input.script} \
                {input.drplt} \
                {input.cvrg} \
                {params.species_id} \
                {params.threshold_horizontal_coverage} \
                {params.threshold_dissimilarity} \
                {input.sample_map} \
                {output.sample_to_scg} \
                {output.scg}
        """


rule filter_metagenotype:
    output:
        "data/{stemA}.metagenotype{stemB}filt-poly{poly}-cvrg{cvrg}.nc",
    input:
        "data/{stemA}.metagenotype{stemB}nc",
    wildcard_constraints:
        poly="[0-9]+",
        cvrg="[0-9]+",
    params:
        poly=lambda w: float(w.poly) / 100,
        cvrg=lambda w: float(w.cvrg) / 100,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts filter_mgen --min-minor-allele-freq {params.poly} --min-horizontal-cvrg {params.cvrg} {input} {output}
        """


rule filter_metagenotype_and_subsample:
    output:
        "data/{stemA}.metagenotype{stemB}filt-poly{poly}-cvrg{cvrg}-g{npos}.nc",
    input:
        "data/{stemA}.metagenotype{stemB}nc",
    wildcard_constraints:
        poly="[0-9]+",
        cvrg="[0-9]+",
        npos="[0-9]+",
    params:
        poly=lambda w: float(w.poly) / 100,
        cvrg=lambda w: float(w.cvrg) / 100,
        npos=lambda w: int(w.npos),
        seed=0,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts filter_mgen --min-minor-allele-freq {params.poly} --min-horizontal-cvrg {params.cvrg} --num-positions {params.npos} --random-seed {params.seed} {input} {output}
        """


checkpoint extract_metagenotype_dimensions:
    output:
        "{stemA}.metagenotype{stemB}.mgen_dims.tsv",
    input:
        script="scripts/extract_metagenotype_dimensions.py",
        mgen="{stemA}.metagenotype{stemB}.nc",
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.mgen} > {output}
        """


localrules:
    extract_metagenotype_dimensions,


# NOTE: This amended function only works if extract_metagenotype_dimensions has already been run.
# This can be setup automatically by making extract_metagenotype_dimensions a checkpoint rule.
def checkpoint_extract_metagenotype_dimensions(wildcards):
    with open(
        checkpoints.extract_metagenotype_dimensions.get(**wildcards).output[0]
    ) as f:
        sizes = {}
        for line in f:
            dim, length = line.split()
            sizes[dim] = int(length)
    return sizes


rule simulate_from_model_no_missing:
    output:
        "data/sfacts_simulate-model_{model_name}-n{n}-g{g}-s{s}-rho{rho_hyper}-pi{pi_hyper}-mu{mu_hyper_mean}-eps{epsilon_hyper_mode}-alpha{alpha_hyper_mean}-seed{seed}.world.nc",
    wildcard_constraints:
        seed="[0-9]+",
        alpha_hyper_mean="[0-9]+",
    params:
        seed=lambda w: int(w.seed),
        n=lambda w: int(w.n),
        g=lambda w: int(w.g),
        s=lambda w: int(w.s),
        rho_hyper=lambda w: float(w.rho_hyper) / 10,
        pi_hyper=lambda w: float(w.pi_hyper) / 100,
        epsilon_hyper_mode=lambda w: float(w.epsilon_hyper_mode) / 1000,
        alpha_hyper_mean=lambda w: float(w.alpha_hyper_mean),
        mu_hyper_mean=lambda w: float(w.mu_hyper_mean) / 10,
        gamma_hyper=1e-5,
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        dd(
            r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output}
        python3 -m sfacts simulate \
                --model-structure {wildcards.model_name} \
                -n {params.n} -g {params.g} -s {params.s} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} pi_hyper={params.pi_hyper} \
                --hyperparameters epsilon_hyper_mode={params.epsilon_hyper_mode} epsilon_hyper_spread=1e3 \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} alpha_hyper_scale=1e-5 \
                --hyperparameters mu_hyper_mean={params.mu_hyper_mean} mu_hyper_scale=1e-5 \
                --hyperparameters m_hyper_r_mean=10. m_hyper_r_scale=1e-5 \
                --seed {params.seed} \
                --outpath {output}
        """
        )


rule simulate_from_simplest_model:
    output:
        "data/sfacts_simulate-model_simplest_simulation-n{n}-g{g}-s{s}-pi{pi_hyper}-mu{mu_hyper_mean}-eps{epsilon_hyper_mode}-seed{seed}.world.nc",
    wildcard_constraints:
        seed="[0-9]+",
        epsilon_hyper_mode="[0-9]+",
    params:
        model_name="simplest_simulation",
        seed=lambda w: int(w.seed),
        n=lambda w: int(w.n),
        g=lambda w: int(w.g),
        s=lambda w: int(w.s),
        pi_hyper=lambda w: float(w.pi_hyper) / 100,
        epsilon_hyper_mode=lambda w: float(w.epsilon_hyper_mode) / 1000,
        mu_hyper_mean=lambda w: float(w.mu_hyper_mean) / 10,
    conda:
        "conda/sfacts.yaml"
    shell:
        dd(
            r"""
        rm -rf {output}
        python3 -m sfacts simulate \
                --model-structure {params.model_name} \
                -n {params.n} -g {params.g} -s {params.s} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters epsilon_hyper_mode={params.epsilon_hyper_mode} \
                --hyperparameters mu_hyper_mean={params.mu_hyper_mean} \
                --seed {params.seed} \
                --outpath {output}
        """
        )


rule extract_and_portion_metagenotype:
    output:
        tsv="{stem}.metagenotype-n{n}-g{g}.tsv",
        nc="{stem}.metagenotype-n{n}-g{g}.nc",
    input:
        script="scripts/extract_metagenotype.py",
        world="{stem}.world.nc",
    wildcard_constraints:
        n="[0-9]+",
        g="[0-9]+",
    params:
        num_samples=lambda w: int(w.n),
        num_positions=lambda w: int(w.g),
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        {input.script} {input.world} {params.num_samples} {params.num_positions} {output.tsv} {output.nc}
        """


localrules:
    extract_and_portion_metagenotype,


rule build_world_from_sfinder_tsv:
    output:
        "{stem}.fit-sfinder{params}.world.nc",
    input:
        script="scripts/sfacts_world_from_flatfiles.py",
        gamma="{stem}.fit-sfinder{params}.gamma.tsv",
        pi="{stem}.fit-sfinder{params}.pi.tsv",
        mgen="{stem}.tsv",
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        {input.script} {input.mgen} {input.gamma} {input.pi} {output}
        """


localrules:
    build_world_from_tsv,


rule evaluate_fit_against_simulation:
    output:
        "data/sfacts_simulate-{sim_stem}.metagenotype-{portion_stem}.fit-{params}.evaluation.tsv",
    input:
        script="scripts/evaluate_haplotyping_against_simulation.py",
        sim="data/sfacts_simulate-{sim_stem}.world.nc",
        fit="data/sfacts_simulate-{sim_stem}.metagenotype-{portion_stem}.fit-{params}.world.nc",
        # bench="data/sfacts_simulate-{sim_stem}.metagenotype-{portion_stem}.fit-{params}.benchmark",
    conda:
        "conda/sfacts.yaml"
    resources:
        walltime_hr=12,
    shell:
        """
        {input.script} {input.sim} {input.fit} {output}
        """

rule compare_inferences_to_scgs:
    output:
        "data/ucfmt.sp-{species_id}.metagenotype.filt-{filt_stem}.fit-{fit_stem}.scg_comparison.tsv",
    input:
        script="scripts/compare_inferred_genotypes_to_scgs.py",
        mgen="data/ucfmt.sp-{species_id}.metagenotype.filt-{filt_stem}.nc",
        scg="data/ucfmt.sp-{species_id}.derep.genotype.nc",
        scg_to_sample="data/ucfmt.sp-{species_id}.derep.barcode_to_sample.tsv",
        library_to_sample="data/ucfmt.barcode_to_sample.tsv",
        fit="data/ucfmt.sp-{species_id}.metagenotype.filt-{filt_stem}.fit-{fit_stem}.world.nc",
    params:
        rabund_threshold=0.01,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.mgen} {input.scg} {input.fit} {input.scg_to_sample} {input.library_to_sample} {params.rabund_threshold} {output}
        """


rule combine_scg_comparisons_for_all_species:
    output:
        "data/ucfmt.filt-{filt_stem}.fit-{fit_stem}.all_scg_comparison.tsv",
    input:
        # lambda w: [
        #     f"data/ucfmt.sp-{s}.metagenotype.filt-{w.filt_stem}.fit-{w.fit_stem}.world.nc"
        #     for s in checkpoint_extract_scg_coverage_table(w, threshold=0.01)
        # ],
        lambda w: [
            f"data/ucfmt.sp-{s}.metagenotype.filt-{w.filt_stem}.fit-{w.fit_stem}.scg_comparison.tsv"
            for s in checkpoint_extract_scg_coverage_table(w, threshold=0.01)
        ],
    params:
        all_species=lambda w: checkpoint_extract_scg_coverage_table(w, threshold=0.01),
        header="\t".join(
            [
                "species_id",
                "scg",
                "sample",
                "mgen",
                "mgen_entropy",
                "mgen_horizontal_coverage",
                "scg_entropy",
                "scg_horizontal_coverage",
                "comm_entropy",
                "fdist_any_strain",
                "ddist_any_strain",
                "mdist_any_strain",
                "fdist_any_mgen",
                "ddist_any_mgen",
                "mdist_any_mgen",
                "fdist_focal_strain",
                "ddist_focal_strain",
                "mdist_focal_strain",
                "fdist_focal_mgen",
                "ddist_focal_mgen",
                "mdist_focal_mgen",
            ]
        ),
        prefix="data/ucfmt.sp-",
        suffix=lambda w: f".metagenotype.filt-{w.filt_stem}.fit-{w.fit_stem}.scg_comparison.tsv",
    shell:
        """
        echo '{params.header}' > {output}
        for s in {params.all_species}
        do
            awk -v species_id=$s -v OFS='\t' 'NR>1{{print species_id,$0}}' {params.prefix}$s{params.suffix}
        done >> {output}
        """


rule recombine_sfacts_genotypes_and_community:
    output:
        "{stemA}.metagenotype{stemB}.fit-{fit_params}.genotype_{genofit_type}-g{nposition}-seed{seed}.world.nc",
    input:
        metagenotype="{stemA}.metagenotype{stemB}.nc",
        community="{stemA}.metagenotype{stemB}.fit-{fit_params}.world_collapsed.nc",
        genotype_chunks=lambda w: [
            f"{{stemA}}.metagenotype{{stemB}}.fit-{{fit_params}}.genotype_{{genofit_type}}-g{{nposition}}-chunk{chunk_i}-seed{{seed}}.nc"
            for chunk_i in range(
                math.ceil(
                    checkpoint_extract_metagenotype_dimensions(w)["position"]
                    / int(w.nposition)
                )
            )
        ],
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output}
        python3 -m sfacts concatenate_genotype_chunks \
                --verbose \
                --community {input.community} --metagenotype {input.metagenotype} \
                --outpath {output} \
                {input.genotype_chunks}
        """


rule evaluate_simulation_fits_at_fixed_s_to_n_ratio:
    output:
        touch(
            "data/evaluate_{fit_type}_at_{n_ratio}x_samples-g{g}-s{s}-{sim_params}-seed{seed}.flag"
        ),
    wildcard_constraints:
        s="[0-9]+",
    input:
        lambda w: [
            "data/sfacts_simulate-model_simplest_simulation-n{n}-g{{g}}-s{{s}}-{{sim_params}}.metagenotype-n{n}-g{{g}}.fit-{{fit_type}}-s{fit_s}-seed{{seed}}.evaluation.tsv".format(
                n=int(w.s) * int(w.n_ratio),
                fit_s=fit_s,
            )
            for fit_s in [
                int(w.s),
                int(w.s) * 2,
            ]
        ],


rule fit_sfacts_strategy44:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts44-s{nstrain}-g{nposition}-seed{seed}.world.nc",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts44-s{nstrain}-g{nposition}-seed{seed}.benchmark"
    wildcard_constraints:
        nstrain="[0-9]+",
        seed="[0-9]+",
        nposition="[0-9]+",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=1e1,
        alpha_hyper_scale=1e-6,
        anneal_gamma_hyper=1e0,
        anneal_rho_hyper=5.0,
        anneal_alpha_hyper_mean=1e1,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=10_000,
        anneal_wait=2000,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.fit}
        python3 -m sfacts fit -m {params.model_name}  \
                --verbose --device {resources.device} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --num-strains {params.nstrain} --num-positions {params.nposition} \
                --nmf-init \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --anneal-hyperparameters gamma_hyper={params.anneal_gamma_hyper} rho_hyper={params.anneal_rho_hyper} pi_hyper=1.0 alpha_hyper_mean={params.anneal_alpha_hyper_mean} \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                {input.data} \
                {output.fit}
        """

rule fit_sfacts_strategy44_ratio_strains:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts44_v-s{rstrain}-g{nposition}-seed{seed}.world.nc",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts44_v-s{rstrain}-g{nposition}-seed{seed}.benchmark"
    wildcard_constraints:
        rstrain="[0-9]+",
        seed="[0-9]+",
        nposition="[0-9]+",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        rstrain=lambda w: float(w.rstrain) / 100,
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=1e1,
        alpha_hyper_scale=1e-6,
        anneal_gamma_hyper=1e0,
        anneal_rho_hyper=5.0,
        anneal_alpha_hyper_mean=1e1,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=10_000,
        anneal_wait=2000,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.fit}
        python3 -m sfacts fit -m {params.model_name}  \
                --verbose --device {resources.device} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --strains-per-sample {params.rstrain} --num-positions {params.nposition} \
                --nmf-init \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --anneal-hyperparameters gamma_hyper={params.anneal_gamma_hyper} rho_hyper={params.anneal_rho_hyper} pi_hyper=1.0 alpha_hyper_mean={params.anneal_alpha_hyper_mean} \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                {input.data} \
                {output.fit}
        """

use rule fit_sfacts_strategy44 as fit_sfacts_strategy44_gpu with:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts44_gpu-s{nstrain}-g{nposition}-seed{seed}.world.nc",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts44_gpu-s{nstrain}-g{nposition}-seed{seed}.benchmark"
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=1e1,
        alpha_hyper_scale=1e-6,
        anneal_gamma_hyper=1e0,
        anneal_rho_hyper=5.0,
        anneal_alpha_hyper_mean=1e1,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=10_000,
        anneal_wait=2000,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device='cuda',
        gpu_mem_mb=5_000,

use rule fit_sfacts_strategy44 as fit_sfacts_strategy44_cpu with:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts44_cpu-s{nstrain}-g{nposition}-seed{seed}.world.nc",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts44_cpu-s{nstrain}-g{nposition}-seed{seed}.benchmark"
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=1e1,
        alpha_hyper_scale=1e-6,
        anneal_gamma_hyper=1e0,
        anneal_rho_hyper=5.0,
        anneal_alpha_hyper_mean=1e1,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=10_000,
        anneal_wait=2000,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        walltime_hr=72,
        pmem=5_000,
        mem_mb=5_000,
        device='cpu',

rule fit_sfacts_strategy44_timeit:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts44_timeit-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        time="{stem}.metagenotype{stemB}.fit-sfacts44_timeit-s{nstrain}-g{nposition}-seed{seed}.time",
    wildcard_constraints:
        nstrain="[0-9]+",
        seed="[0-9]+",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=1e1,
        alpha_hyper_scale=1e-6,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        walltime_hr=1,
        pmem=5_000,
        mem_mb=5_000,
        device="cpu",
        gpu_mem_mb=0,
    conda:
        "conda/sfacts.yaml"
    shell:
        r"""
        rm -rf {output.fit}
        `which time` -o {output.time} -- \
            python3 -m sfacts fit -m {params.model_name}  \
                --verbose --device {resources.device} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --num-strains {params.nstrain} --num-positions {params.nposition} \
                --nmf-init \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 10 --lag1 {params.lag1} --lag2 {params.lag2} \
                {input.data} \
                {output.fit}
        """


rule fit_sfacts_strategy44_gpumem:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts44_gpumem-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        gpumem="{stem}.metagenotype{stemB}.fit-sfacts44_gpumem-s{nstrain}-g{nposition}-seed{seed}.gpumem",
    wildcard_constraints:
        nstrain="[0-9]+",
        seed="[0-9]+",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=1e1,
        alpha_hyper_scale=1e-6,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        lr=0.05,
        min_learning_rate=1e-2,
    resources:
        pmem=5_000,
        mem_mb=5_000,
        device="cuda",
        gpu_mem_mb=3_000,
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.fit}
        nvidia-smi -i $CUDA_VISIBLE_DEVICES --query-gpu=memory.used --format=csv,noheader,nounits --loop-ms=100 --filename={output.gpumem} &
        gpumem_pid=$!
            python3 -m sfacts fit -m {params.model_name}  \
                --verbose --device {resources.device} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --num-strains {params.nstrain} --num-positions {params.nposition} \
                --nmf-init \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 1_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                {input.data} \
                {output.fit}
        sleep 5
        kill $gpumem_pid
        wait -n $gpumem_pid
        """


# NOTE: The below 'drop_*' rules are included to bring sfacts fits into
# alignment with sfinder fit file naming, while keeping the
# npositions-filename-parameterization which is needed for other workflows.
rule drop_g1000000_from_world_suffix:
    output: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-seed{seed}.world.nc'
    input: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-g1000000-seed{seed}.world.nc'
    wildcard_constraints:
        nstrain='[0-9]+',
        seed='[0-9]+',
    shell: alias_recipe


localrules: drop_g1000000_from_world_suffix


rule drop_g1000000_from_benchmark_suffix:
    output: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-seed{seed}.benchmark'
    input: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-g1000000-seed{seed}.benchmark'
    shell: alias_recipe


localrules: drop_g1000000_from_benchmark_suffix


rule drop_g1000000_from_time_suffix:
    output: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-seed{seed}.time'
    input: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-g1000000-seed{seed}.time'
    wildcard_constraints:
        nstrain='[0-9]+',
        seed='[0-9]+',
    shell: alias_recipe


localrules: drop_g1000000_from_time_suffix


rule drop_g1000000_from_gpumem_suffix:
    output: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-seed{seed}.gpumem'
    input: '{stem}.fit-sfacts{_fit_type}-s{nstrain}-g1000000-seed{seed}.gpumem'
    wildcard_constraints:
        nstrain='[0-9]+',
        seed='[0-9]+',
    shell: alias_recipe


localrules: drop_g1000000_from_gpumem_suffix


rule fit_sfacts_strategy41_genotypes:
    output:
        genotype_chunk_fit="{stem}.metagenotype{stemB}.fit-{fit_params}.refit-sfacts41-g{nposition}-chunk{chunk_i}-seed{seed}.nc",
    input:
        metagenotype="{stem}.metagenotype{stemB}.nc",
        community="{stem}.metagenotype{stemB}.fit-{fit_params}.world.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nposition=lambda w: int(w.nposition),
        npositionB=lambda w: min(int(w.nposition), max(int(w.nposition) // 5, 500)),
        block_number=lambda w: int(w.chunk_i),
        precision=64,
        gamma_hyper=1.0,
        alpha_hyper_mean=200,
        alpha_hyper_scale=1e-6,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        lr=0.5,
    resources:
        # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
        # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
        walltime_hr=12,
        pmem=5_000,
        mem_mb=5_000,
        gpu_mem_mb={0: 0, 1: 4_000}[config["USE_CUDA"]],
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.genotype_chunk_fit}
        python3 -m sfacts fit_genotype -m {params.model_name}  \
                --num-positions {params.nposition} --block-number {params.block_number} --num-positionsB {params.npositionB} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --optimizer-learning-rate {params.lr} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --verbose --device {params.device} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                {input.community} {input.metagenotype} {output.genotype_chunk_fit}
        """


rule recombine_refit_genotypes_and_community:
    output:
        "{stemA}.metagenotype{stemB}.fit-{fit_params}.refit-{refit_type}-g{nposition}-seed{seed}.world.nc",
    input:
        metagenotype="{stemA}.metagenotype{stemB}.nc",
        community="{stemA}.metagenotype{stemB}.fit-{fit_params}.world.nc",
        genotype_chunks=lambda w: [
            f"{{stemA}}.metagenotype{{stemB}}.fit-{{fit_params}}.refit-{{refit_type}}-g{{nposition}}-chunk{chunk_i}-seed{{seed}}.nc"
            for chunk_i in range(
                math.ceil(
                    checkpoint_extract_metagenotype_dimensions(w)["position"]
                    / int(w.nposition)
                )
            )
        ],
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        rm -rf {output}
        python3 -m sfacts concatenate_genotype_chunks \
                --verbose \
                --community {input.community} --metagenotype {input.metagenotype} \
                --outpath {output} \
                {input.genotype_chunks}
        """
