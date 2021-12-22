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
    output: "{stemA}.metagenotype{stemB}.mgen_dims.tsv"
    input:
        script="scripts/extract_metagenotype_dimensions.py",
        mgen="{stemA}.metagenotype{stemB}.nc",
    conda: "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.mgen} > {output}
        """


localrules:
    extract_metagenotype_dimensions,


def checkpoint_extract_metagenotype_dimensions(wildcards):
    with open(checkpoints.extract_metagenotype_dimensions.get(**wildcards).output[0]) as f:
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
        model_name='simplest_simulation',
        seed=lambda w: int(w.seed),
        n=lambda w: int(w.n),
        g=lambda w: int(w.g),
        s=lambda w: int(w.s),
        pi_hyper=lambda w: float(w.pi_hyper) / 100,
        epsilon_hyper_mode=lambda w: float(w.epsilon_hyper_mode) / 1000,
        mu_hyper_mean=lambda w: float(w.mu_hyper_mean) / 10,
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        dd(
            r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
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


# rule fit_sfacts_strategy1:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts1-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts1-s{nstrain}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype-n{n}-g{g}.tsv",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=0.0005,
#         rho_hyper=0.08,
#         pi_hyper=0.2,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=500,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1-s{nstrain}-seed{seed}.benchmark"
#     resources:
#         pmem=resource_calculator(data=20, nstrain=1, agg=math.prod),
#         gpu_mem_mb=resource_calculator(data=20, nstrain=1, agg=math.prod),
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         """
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_simple -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --inpath {input.data} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --optimizer {params.optimizer} --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 -s {params.nstrain} \
#                 --random-seed {params.seed} \
#                 --outpath {output.world} --history-out {output.history}
#         """
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy2 with:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts2-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts2-s{nstrain}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=1e-6,
#         rho_hyper=1e5,
#         pi_hyper=1.0,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=100,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts2-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy4 with:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts4-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts4-s{nstrain}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=1e-6,
#         rho_hyper=1e5,
#         pi_hyper=1.0,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=100,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts4-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy3 with:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts3-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts3-s{nstrain}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=1e-6,
#         rho_hyper=1e5,
#         pi_hyper=1.0,
#         seed=lambda w: int(w.seed),
#         model_name="simple",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=100,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts3-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy5 with:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts5-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts5-s{nstrain}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=0.0005,
#         rho_hyper=0.08,
#         pi_hyper=0.2,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=100,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts5-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy6 with:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts6-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts6-s{nstrain}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=0.00005,
#         rho_hyper=0.08,
#         pi_hyper=0.3,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=500,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts6-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy7 with:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts7-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts7-s{nstrain}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=0.5,
#         rho_hyper=1.0,
#         pi_hyper=0.5,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=500,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts7-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategyN with:
#     output:
#         world="{stem}.metagenotype-n{n}-g{g}.fit-sfactsN-pi{pi_hyper}-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype-n{n}-g{g}.fit-sfactsN-pi{pi_hyper}-s{nstrain}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=0.0005,
#         rho_hyper=0.08,
#         pi_hyper=lambda w: float(w.pi_hyper) / 100,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         optimizer="Adamax",
#         lag1=50,
#         lag2=500,
#         lr=5e-1,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfactsN-pi{pi_hyper}-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy1_cpu with:
#     output:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_cpu-s{nstrain}-seed{seed}.world.nc",
#     params:
#         device="cpu",
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=0.0005,
#         rho_hyper=0.08,
#         pi_hyper=0.2,
#         learning_rate=1e-2,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2",
#         optimizer="Adamax",
#         lag1=40,
#         lag2=200,
#         lr=1e-3,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_cpu-s{nstrain}-seed{seed}.benchmark"
#
#
# use rule fit_sfacts_strategy1 as fit_sfacts_strategy1_gpu with:
#     output:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_gpu-s{nstrain}-seed{seed}.world.nc",
#     params:
#         device="cuda",
#         nstrain=lambda w: int(w.nstrain),
#         gamma_hyper=0.0005,
#         rho_hyper=0.08,
#         pi_hyper=0.2,
#         learning_rate=1e-2,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2",
#         optimizer="Adamax",
#         lag1=40,
#         lag2=200,
#         lr=1e-3,
#     benchmark:
#         "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_gpu-s{nstrain}-seed{seed}.benchmark"


rule build_world_from_sfinder_tsv:
    output:
        "{stem}.fit-{params}.world.nc",
    input:
        script="scripts/sfacts_world_from_flatfiles.py",
        gamma="{stem}.fit-{params}.gamma.tsv",
        pi="{stem}.fit-{params}.pi.tsv",
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
        bench="data/sfacts_simulate-{sim_stem}.metagenotype-{portion_stem}.fit-{params}.benchmark",
    # conda:
    #     "conda/sfacts.yaml"
    resources:
        walltime_hr=12,
    shell:
        """
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        {input.script} {input.sim} {input.fit} {input.bench} {output}
        """
#
#
# rule fit_sfacts_strategy8:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts8-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts8-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         precision=32,
#         gamma_hyper=1e-5,
#         rho_hyper=0.1,
#         pi_hyper=0.2,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2",
#         lag1=50,
#         lag2=100,
#         anneal_steps=2500,
#         anneal_wait=500,
#         lr=0.1,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         pmem=resource_calculator(data=2, nposition=0.002, nstrain=2, agg=math.prod),
#         gpu_mem_mb=resource_calculator(
#             data=2, nposition=0.002, nstrain=2, agg=math.prod
#         ),
#         walltime_hr=12,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_complex -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --anneal-hyperparameters gamma_hyper=1.0 pi_hyper=1.0 rho_hyper=1.0 \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --refinement-hyperparameters gamma_hyper=1.0 rho_hyper=5.0 pi_hyper=0.1 \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --precision {params.precision} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
#
# use rule fit_sfacts_strategy8 as fit_sfacts_strategy9 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts9-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts9-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         precision=32,
#         gamma_hyper=1e-5,
#         rho_hyper=0.1,
#         pi_hyper=0.2,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2",
#         lag1=50,
#         lag2=500,
#         anneal_steps=5000,
#         anneal_wait=500,
#         lr=0.1,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         pmem=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
#         gpu_mem_mb=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
#         walltime_hr=24,
#
#
# use rule fit_sfacts_strategy8 as fit_sfacts_strategy10 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts10-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts10-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         precision=32,
#         gamma_hyper=1e-3,
#         rho_hyper=0.1,
#         pi_hyper=0.2,
#         seed=lambda w: int(w.seed),
#         model_name="simple",
#         lag1=50,
#         lag2=500,
#         anneal_steps=5000,
#         anneal_wait=500,
#         lr=0.1,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         pmem=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
#         gpu_mem_mb=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
#         walltime_hr=24,
#
#
# rule fit_sfacts_strategy11:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts11-s{nstrain}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts11-s{nstrain}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         precision=32,
#         gamma_hyper=1e-4,
#         rho_hyper=0.1,
#         pi_hyper=0.2,
#         alpha_hyper_scale=0.1,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=2500,
#         anneal_wait=500,
#         lr=0.05,
#     resources:
#         pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_simple -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --anneal-hyperparameters gamma_hyper=1.0 pi_hyper=1.0 rho_hyper=1.0 \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
#
# rule fit_sfacts_strategy12:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts12-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts12-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         precision=32,
#         gamma_hyper=1e-4,
#         rho_hyper=0.1,
#         pi_hyper=0.2,
#         alpha_hyper_scale=0.1,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=2500,
#         anneal_wait=500,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_complex -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --refinement-hyperparameters gamma_hyper=1.0 \
#                 --anneal-hyperparameters gamma_hyper=1.0 pi_hyper=1.0 rho_hyper=1.0 \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
#
# use rule fit_sfacts_strategy12 as fit_sfacts_strategy13 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts13-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts13-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         precision=64,
#         gamma_hyper=1e-2,
#         rho_hyper=1.0,
#         pi_hyper=0.2,
#         alpha_hyper_scale=0.1,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=2500,
#         anneal_wait=500,
#         lr=0.001,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         pmem=5000,
#         gpu_mem_mb=5000,
#         walltime_hr=48,
#
#
# use rule fit_sfacts_strategy12 as fit_sfacts_strategy14 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts14-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts14-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         precision=32,
#         gamma_hyper=5e-5,
#         rho_hyper=0.1,
#         pi_hyper=0.2,
#         alpha_hyper_scale=0.1,
#         seed=lambda w: int(w.seed),
#         model_name="simple_ssdd2_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=2500,
#         anneal_wait=500,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
#
# rule fit_sfacts_strategy15:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts15-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts15-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-10,
#         rho_hyper=0.3,
#         pi_hyper=0.5,
#         alpha_hyper_mean=200,
#         alpha_hyper_scale=0.1,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_complex -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --refinement-hyperparameters gamma_hyper=1.0 \
#                 --anneal-hyperparameters pi_hyper=1.0 rho_hyper=1.0 \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --num-positionsB {params.npositionB} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
# use rule fit_sfacts_strategy15 as fit_sfacts_strategy16 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts16-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts16-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-10,
#         rho_hyper=0.3,
#         pi_hyper=0.5,
#         alpha_hyper_mean=2e2,
#         alpha_hyper_scale=2.0,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,


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
        threshold=0.01,
        pseudo=1e-10,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.mgen} {input.scg} {input.fit} {input.scg_to_sample} {input.library_to_sample} {params.threshold} {params.pseudo} {output}
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
                "focal_sample",
                "smallest_fdist_all_strains",
                "smallest_ddist_all_strains",
                "smallest_fdist_all_mgen",
                "smallest_ddist_all_mgen",
                "smallest_fdist_focal_strains",
                "smallest_ddist_focal_strains",
                "smallest_fdist_focal_mgen",
                "smallest_ddist_focal_mgen",
                "scg_horizontal_coverage",
                "horizontal_coverage_focal_mgen",
                "community_entropy_focal_sample",
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


# localrules:
#     run_scg_comparison_all_species,
#
#
# use rule fit_sfacts_strategy15 as fit_sfacts_strategy17 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts17-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts17-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition),
#         precision=64,
#         gamma_hyper=1e-10,
#         rho_hyper=1.0,
#         pi_hyper=0.2,
#         alpha_hyper_mean=2e2,
#         alpha_hyper_scale=2.0,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.001,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=168,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#
# rule fit_sfacts_strategy18:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts18-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts18-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-10,
#         rho_hyper=1.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=200,
#         alpha_hyper_scale=0.1,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_complex -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --refinement-hyperparameters gamma_hyper=1.0 \
#                 --anneal-hyperparameters pi_hyper=1.0 rho_hyper=1.0 \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --num-positionsB {params.npositionB} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
# use rule fit_sfacts_strategy18 as fit_sfacts_strategy19 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts19-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts19-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-10,
#         rho_hyper=1.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=10,
#         alpha_hyper_scale=0.01,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#
# use rule fit_sfacts_strategy18 as fit_sfacts_strategy20 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts20-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts20-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-15,
#         rho_hyper=10.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=2e2,
#         alpha_hyper_scale=2.0,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#
# use rule fit_sfacts_strategy18 as fit_sfacts_strategy21 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts21-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts21-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-15,
#         rho_hyper=10.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=1e6,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#
# rule fit_sfacts_strategy22:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts22-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts22-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-10,
#         rho_hyper=1.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=2000,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_complex -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --refinement-hyperparameters gamma_hyper=1.0 \
#                 --anneal-hyperparameters pi_hyper=1.0 rho_hyper=1.0 alpha_hyper_mean=1.0 \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --num-positionsB {params.npositionB} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
#
# use rule fit_sfacts_strategy22 as fit_sfacts_strategy23 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts23-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts23-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-15,
#         rho_hyper=0.5,
#         pi_hyper=0.3,
#         alpha_hyper_mean=100,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
#
# use rule fit_sfacts_strategy22 as fit_sfacts_strategy24 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts24-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts24-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=1.0,
#         pi_hyper=0.3,
#         alpha_hyper_mean=200,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# use rule fit_sfacts_strategy22 as fit_sfacts_strategy25 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts25-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts25-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=1.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=200,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# use rule fit_sfacts_strategy22 as fit_sfacts_strategy26 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts26-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts26-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=10.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=200,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# use rule fit_sfacts_strategy22 as fit_sfacts_strategy27 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts27-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts27-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=1.0,
#         pi_hyper=0.3,
#         alpha_hyper_mean=100,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# use rule fit_sfacts_strategy22 as fit_sfacts_strategy28 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts28-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts28-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=2.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=1000,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# use rule fit_sfacts_strategy22 as fit_sfacts_strategy29 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts29-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts29-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=2.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=100,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# rule fit_sfacts_strategy30:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts30-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts30-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=2.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=1e3,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_complex -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --refinement-hyperparameters gamma_hyper=1.0 \
#                 --nmf-init \
#                 --anneal-hyperparameters pi_hyper=1.0 alpha_hyper_mean=1.0 \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --num-positionsB {params.npositionB} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
# use rule fit_sfacts_strategy30 as fit_sfacts_strategy31 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts31-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts31-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=2.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=1e3,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# use rule fit_sfacts_strategy30 as fit_sfacts_strategy32 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts32-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts32-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-3,
#         rho_hyper=5.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=1e3,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=9000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# use rule fit_sfacts_strategy30 as fit_sfacts_strategy33 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts33-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts33-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-3,
#         rho_hyper=0.5,
#         pi_hyper=0.3,
#         alpha_hyper_mean=2e2,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=9000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#
# rule fit_sfacts_strategy34:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts34-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts34-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-10,
#         rho_hyper=0.5,
#         anneal_rho_hyper=5.0,
#         pi_hyper=0.3,
#         alpha_hyper_mean=1e3,
#         anneal_alpha_hyper_mean=1.0,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=9000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.history} {output.world}
#         python3 -m sfacts fit_complex -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --refinement-hyperparameters gamma_hyper=1.0 \
#                 --nmf-init \
#                 --anneal-hyperparameters rho_hyper={params.anneal_rho_hyper} pi_hyper=1.0 alpha_hyper_mean={params.anneal_alpha_hyper_mean} \
#                 --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --num-positionsB {params.npositionB} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --history-out {output.history} \
#                 {input.data} {output.world}
#         """
#
# use rule fit_sfacts_strategy34 as fit_sfacts_strategy35 with:
#     output:
#         world="{stem}.metagenotype{stemB}.fit-sfacts35-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         history="{stem}.metagenotype{stemB}.fit-sfacts35-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-10,
#         rho_hyper=0.5,
#         pi_hyper=0.3,
#         alpha_hyper_mean=10.,
#         alpha_hyper_scale=1e-6,
#         anneal_rho_hyper=5.0,
#         anneal_alpha_hyper_mean=1e3,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=9000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,

# NOTE: This is my current favorite for matching SCGs
rule fit_sfacts_strategy36:
    output:
        initial_fit="{stem}.metagenotype{stemB}.fit-sfacts36-s{nstrain}-g{nposition}-seed{seed}.world_initial.nc",
        collapsed_fit="{stem}.metagenotype{stemB}.fit-sfacts36-s{nstrain}-g{nposition}-seed{seed}.world_collapsed.nc",
        full_fit="{stem}.metagenotype{stemB}.fit-sfacts36-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        # history="{stem}.metagenotype{stemB}.fit-sfacts36-s{nstrain}-g{nposition}-seed{seed}.history.list",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts36-s{nstrain}-g{nposition}-seed{seed}.benchmark",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        npositionB=lambda w: int(w.nposition) * 10,
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=10.,
        alpha_hyper_scale=1e-6,
        anneal_rho_hyper=5.0,
        anneal_alpha_hyper_mean=1e3,
        refine_alpha_hyper_mean=200,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=4000,
        anneal_wait=1000,
        lr=0.05,
        collapse=0.02,
        cull=0.001,
    resources:
        # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
        # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.initial_fit} {output.collapsed_fit} {output.full_fit}
        python3 -m sfacts fit_complex2 -m {params.model_name}  \
                --verbose --device {resources.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --nmf-init \
                --refinement-hyperparameters gamma_hyper=1.0 alpha_hyper_mean={params.refine_alpha_hyper_mean} \
                --anneal-hyperparameters rho_hyper={params.anneal_rho_hyper} pi_hyper=1.0 alpha_hyper_mean={params.anneal_alpha_hyper_mean} \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --optimizer-learning-rate {params.lr} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --num-positions {params.nposition} \
                --num-positionsB {params.npositionB} \
                --collapse {params.collapse} --cull {params.cull} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --outpath0 {output.initial_fit} --outpath1 {output.collapsed_fit} --outpath2 {output.full_fit} \
                {input.data}
        """

use rule fit_sfacts_strategy36 as fit_sfacts_strategy36_cpu with:
    output:
        initial_fit="{stem}.metagenotype{stemB}.fit-sfacts36_cpu-s{nstrain}-g{nposition}-seed{seed}.world_initial.nc",
        collapsed_fit="{stem}.metagenotype{stemB}.fit-sfacts36_cpu-s{nstrain}-g{nposition}-seed{seed}.world_collapsed.nc",
        full_fit="{stem}.metagenotype{stemB}.fit-sfacts36_cpu-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        # history="{stem}.metagenotype{stemB}.fit-sfacts36_cpu-s{nstrain}-g{nposition}-seed{seed}.history.list",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts36_cpu-s{nstrain}-g{nposition}-seed{seed}.benchmark",
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device="cpu",

use rule fit_sfacts_strategy36 as fit_sfacts_strategy36_gpu with:
    output:
        initial_fit="{stem}.metagenotype{stemB}.fit-sfacts36_gpu-s{nstrain}-g{nposition}-seed{seed}.world_initial.nc",
        collapsed_fit="{stem}.metagenotype{stemB}.fit-sfacts36_gpu-s{nstrain}-g{nposition}-seed{seed}.world_collapsed.nc",
        full_fit="{stem}.metagenotype{stemB}.fit-sfacts36_gpu-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        # history="{stem}.metagenotype{stemB}.fit-sfacts36_gpu-s{nstrain}-g{nposition}-seed{seed}.history.list",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts36_gpu-s{nstrain}-g{nposition}-seed{seed}.benchmark",
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device="cuda",
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],

# use rule fit_sfacts_strategy36 as fit_sfacts_strategy37 with:
#     output:
#         initial_fit="{stem}.metagenotype{stemB}.fit-sfacts37-s{nstrain}-g{nposition}-seed{seed}.world_initial.nc",
#         collapsed_fit="{stem}.metagenotype{stemB}.fit-sfacts37-s{nstrain}-g{nposition}-seed{seed}.world_collapsed.nc",
#         full_fit="{stem}.metagenotype{stemB}.fit-sfacts37-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         # history="{stem}.metagenotype{stemB}.fit-sfacts37-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=0.5,
#         pi_hyper=0.3,
#         alpha_hyper_mean=100,
#         alpha_hyper_scale=1e-6,
#         anneal_rho_hyper=5.0,
#         anneal_alpha_hyper_mean=1e3,
#         refine_alpha_hyper_mean=100,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,

# rule fit_sfacts_strategy38:
#     output:
#         initial_fit="{stem}.metagenotype{stemB}.fit-sfacts38-s{nstrain}-g{nposition}-seed{seed}.world_initial.nc",
#         collapsed_fit="{stem}.metagenotype{stemB}.fit-sfacts38-s{nstrain}-g{nposition}-seed{seed}.world_collapsed.nc",
#         full_fit="{stem}.metagenotype{stemB}.fit-sfacts38-s{nstrain}-g{nposition}-seed{seed}.world.nc",
#         # history="{stem}.metagenotype{stemB}.fit-sfacts38-s{nstrain}-g{nposition}-seed{seed}.history.list",
#     input:
#         data="{stem}.metagenotype{stemB}.nc",
#     params:
#         device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
#         nstrain=lambda w: int(w.nstrain),
#         nposition=lambda w: int(w.nposition),
#         npositionB=lambda w: int(w.nposition) * 10,
#         precision=32,
#         gamma_hyper=1e-20,
#         rho_hyper=2.0,
#         pi_hyper=0.5,
#         alpha_hyper_mean=1000,
#         alpha_hyper_scale=1e-6,
#         seed=lambda w: int(w.seed),
#         model_name="ssdd3_with_error",
#         lag1=50,
#         lag2=100,
#         anneal_steps=4000,
#         anneal_wait=1000,
#         lr=0.05,
#         collapse=0.02,
#         cull=0.001,
#     resources:
#         # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
#         walltime_hr=12,
#         pmem=5_000,
#         mem_mb=5_000,
#         gpu_mem_mb=5_000,
#     # conda:
#     #     "conda/sfacts.yaml"
#     shell:
#         r"""
#         export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
#         rm -rf {output.initial_fit} {output.collapsed_fit} {output.full_fit}
#         python3 -m sfacts fit_complex2 -m {params.model_name}  \
#                 --verbose --device {params.device} \
#                 --hyperparameters gamma_hyper={params.gamma_hyper} \
#                 --hyperparameters pi_hyper={params.pi_hyper} \
#                 --hyperparameters rho_hyper={params.rho_hyper} \
#                 --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
#                 --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
#                 --nmf-init \
#                 --refinement-hyperparameters gamma_hyper=1.0 \
#                 --optimizer-learning-rate {params.lr} \
#                 --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
#                 --num-strains {params.nstrain} \
#                 --num-positions {params.nposition} \
#                 --num-positionsB {params.npositionB} \
#                 --collapse {params.collapse} --cull {params.cull} \
#                 --precision {params.precision} \
#                 --random-seed {params.seed} \
#                 --outpath0 {output.initial_fit} --outpath1 {output.collapsed_fit} --outpath2 {output.full_fit} \
#                 {input.data}
#         """

# NOTE: This is my current favorite for BIG data.
rule fit_sfacts_strategy39_communities:
    output:
        initial_fit="{stem}.metagenotype{stemB}.fit-sfacts39-s{nstrain}-g{nposition}-seed{seed}.world_initial.nc",
        collapsed_fit="{stem}.metagenotype{stemB}.fit-sfacts39-s{nstrain}-g{nposition}-seed{seed}.world_collapsed.nc",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=64,
        gamma_hyper=1e-10,
        rho_hyper=1.0,
        pi_hyper=0.5,
        alpha_hyper_mean=10.,
        alpha_hyper_scale=1e-6,
        anneal_alpha_hyper_mean=1e3,
        refine_alpha_hyper_mean=200,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=4000,
        anneal_wait=1000,
        lr=0.01,
        collapse=0.02,
        cull=0.001,
    resources:
        # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
        # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
        walltime_hr=168,
        pmem=5_000,
        mem_mb=5_000,
        gpu_mem_mb={0: 0, 1: 10_000}[config["USE_CUDA"]]
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.initial_fit} {output.collapsed_fit}
        python3 -m sfacts fit_community -m {params.model_name}  \
                --verbose --device {params.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --refinement-hyperparameters gamma_hyper=1.0 alpha_hyper_mean={params.refine_alpha_hyper_mean} \
                --anneal-hyperparameters pi_hyper=1.0 alpha_hyper_mean={params.anneal_alpha_hyper_mean} \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --optimizer-learning-rate {params.lr} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --num-positions {params.nposition} \
                --collapse {params.collapse} --cull {params.cull} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --outpath-initial {output.initial_fit} \
                {input.data} {output.collapsed_fit}
        """

rule fit_sfacts_strategy39_genotypes:
    output:
        genotype_chunk_fit="{stem}.metagenotype{stemB}.fit-{fit_params}.genotype_refit39-g{nposition}-chunk{chunk_i}-seed{seed}.nc",
    input:
        metagenotype="{stem}.metagenotype{stemB}.nc",
        community="{stem}.metagenotype{stemB}.fit-{fit_params}.world_collapsed.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nposition=lambda w: int(w.nposition),
        npositionB=lambda w: min(int(w.nposition), max(int(w.nposition) // 10, 500)),
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
        gpu_mem_mb={0: 0, 1: 4_000}[config["USE_CUDA"]]
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

rule recombine_sfacts_genotypes_and_community:
    output:
        "{stemA}.metagenotype{stemB}.fit-{fit_params}.genotype_{genofit_type}-g{nposition}-seed{seed}.world.nc",
    input:
        metagenotype="{stemA}.metagenotype{stemB}.nc",
        community="{stemA}.metagenotype{stemB}.fit-{fit_params}.world_collapsed.nc",
        genotype_chunks=lambda w: [
            f"{{stemA}}.metagenotype{{stemB}}.fit-{{fit_params}}.genotype_{{genofit_type}}-g{{nposition}}-chunk{chunk_i}-seed{{seed}}.nc"
            for chunk_i in range(math.ceil(checkpoint_extract_metagenotype_dimensions(w)['position'] / int(w.nposition)))
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


# NOTE: Same hyperparameters as strategy36 but intended to reflect the speed of the core algorithm.
# No subsampling, only the initial fitting,
# TODO: shorter annealing phase?
rule fit_sfacts_strategy40:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts40-s{nstrain}-seed{seed}.world.nc",
        # history="{stem}.metagenotype{stemB}.fit-sfacts40-s{nstrain}-seed{seed}.history.list",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts40-s{nstrain}-seed{seed}.benchmark",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=10.,
        alpha_hyper_scale=1e-6,
        anneal_rho_hyper=5.0,
        anneal_alpha_hyper_mean=1e3,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=4000,
        anneal_wait=1000,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
        # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.fit}
        python3 -m sfacts fit_community0 -m {params.model_name}  \
                --verbose --device {resources.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --nmf-init \
                --anneal-hyperparameters rho_hyper={params.anneal_rho_hyper} pi_hyper=1.0 alpha_hyper_mean={params.anneal_alpha_hyper_mean} \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                {input.data} \
                {output.fit}
        """

use rule fit_sfacts_strategy40 as fit_sfacts_strategy40_cpu with:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts40_cpu-s{nstrain}-seed{seed}.world.nc",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts40_cpu-s{nstrain}-seed{seed}.benchmark",
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device="cpu",

use rule fit_sfacts_strategy40 as fit_sfacts_strategy40_gpu with:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts40_gpu-s{nstrain}-seed{seed}.world.nc",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts40_gpu-s{nstrain}-seed{seed}.benchmark",
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device='cuda',
        gpu_mem_mb=5_000,

use rule fit_sfacts_strategy40 as fit_sfacts_strategy40_big with:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts40_big-s{nstrain}-seed{seed}.world.nc",
        # history="{stem}.metagenotype{stemB}.fit-sfacts40_big-s{nstrain}-seed{seed}.history.list",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts40_big-s{nstrain}-seed{seed}.benchmark",
    params:
        nstrain=lambda w: int(w.nstrain),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=10.,
        alpha_hyper_scale=1e-6,
        anneal_rho_hyper=5.0,
        anneal_alpha_hyper_mean=1e3,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=4000,
        anneal_wait=1000,
        lr=0.005,
        min_learning_rate=1e-8,


rule evaluate_simulation_fits_at_fixed_s_to_n_ratio:
    output: touch("data/evaluate_{fit_type}_at_{n_ratio}x_samples-g{g}-s{s}-{sim_params}-seed{seed}.flag")
    wildcard_constraints:
        s='[0-9]+'
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
        ]

# RAM profile sfacts40
# FIXME: Rename timeit to memprof
rule fit_sfacts_strategy40_timeit:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts40_timeit-s{nstrain}-seed{seed}.world.nc",
        time="{stem}.metagenotype{stemB}.fit-sfacts40_timeit-s{nstrain}-seed{seed}.time",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts40_timeit-s{nstrain}-seed{seed}.benchmark",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=10.,
        alpha_hyper_scale=1e-6,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        # pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
        # gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
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
            python3 -m sfacts fit_community0 -m {params.model_name}  \
                --verbose --device {resources.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 1_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                {input.data} \
                {output.fit}
        """

rule fit_sfacts_strategy40_gpumem:
    output:
        fit="{stem}.metagenotype{stemB}.fit-sfacts40_gpumem-s{nstrain}-seed{seed}.world.nc",
        gpumem="{stem}.metagenotype{stemB}.fit-sfacts40_gpumem-s{nstrain}-seed{seed}.gpumem",
    benchmark:
        "{stem}.metagenotype{stemB}.fit-sfacts40_gpumem-s{nstrain}-seed{seed}.benchmark",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        precision=32,
        gamma_hyper=1e-10,
        rho_hyper=0.5,
        pi_hyper=0.3,
        alpha_hyper_mean=10.,
        alpha_hyper_scale=1e-6,
        seed=lambda w: int(w.seed),
        model_name="ssdd3_with_error",
        lag1=50,
        lag2=100,
        lr=0.05,
        min_learning_rate=1e-6,
    resources:
        walltime_min=10,
        pmem=5_000,
        mem_mb=5_000,
        device="cuda",
        gpu_mem_mb=10_000,
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.fit}
        nvidia-smi -i $CUDA_VISIBLE_DEVICES --query-gpu=memory.used --format=csv,noheader,nounits --loop-ms=1000 --filename={output.gpumem} &
        gpumem_pid=$!
            python3 -m sfacts fit_community0 -m {params.model_name}  \
                --verbose --device {resources.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_mean={params.alpha_hyper_mean} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 1_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                {input.data} \
                {output.fit}
        sleep 5
        kill $gpumem_pid
        wait -n $gpumem_pid
        """
