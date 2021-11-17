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
    output: 'data/gtprodb.sp-{species}.genotype.nc'
    input:
        script='scripts/vcf_to_gtpro_pileup.py',
        gtpro_snp_dict='ref/gtpro/variants_main.covered.hq.snp_dict.tsv',
        vcf='raw/gtpro_refs/variation_in_species/{species}/core_snps.vcf.gz',
    shell:
        """
        {input.script} {input.gtpro_snp_dict} {input.vcf} {wildcards.species} {output}
        """

rule filter_metagenotype:
    output: "data/{stemA}.metagenotype{stemB}filt-poly{poly}-cvrg{cvrg}.nc"
    input: "data/{stemA}.metagenotype{stemB}nc"
    wildcard_constraints:
        poly='[0-9]+',
        cvrg='[0-9]+',
    params:
        poly=lambda w: float(w.poly) / 100,
        cvrg=lambda w: float(w.cvrg) / 100,
    conda: "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts filter_mgen --min-minor-allele-freq {params.poly} --min-horizontal-cvrg {params.cvrg} {input} {output}
        """

rule filter_metagenotype_and_subsample:
    output: "data/{stemA}.metagenotype{stemB}filt-poly{poly}-cvrg{cvrg}-g{npos}.nc"
    input: "data/{stemA}.metagenotype{stemB}nc"
    wildcard_constraints:
        poly='[0-9]+',
        cvrg='[0-9]+',
        npos='[0-9]+',
    params:
        poly=lambda w: float(w.poly) / 100,
        cvrg=lambda w: float(w.cvrg) / 100,
        npos=lambda w: int(w.npos),
        seed=0,
    conda: "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts filter_mgen --min-minor-allele-freq {params.poly} --min-horizontal-cvrg {params.cvrg} --num-positions {params.npos} --random-seed {params.seed} {input} {output}
        """



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


rule extract_and_portion_metagenotype_tsv:
    output:
        "{stem}.metagenotype-n{n}-g{g}.tsv",
    wildcard_constraints:
        n="[0-9]+",
        g="[0-9]+",
    input:
        script="scripts/extract_metagenotype_to_tsv.py",
        world="{stem}.world.nc",
    params:
        num_samples=lambda w: int(w.n),
        num_positions=lambda w: int(w.g),
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        {input.script} {input.world} {params.num_samples} {params.num_positions} > {output}
        """


localrules:
    extract_and_portion_metagenotype_tsv,


rule fit_sfacts_strategy1:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts1-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts1-s{nstrain}-seed{seed}.history.list",
    input:
        data="{stem}.metagenotype-n{n}-g{g}.tsv",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=0.0005,
        rho_hyper=0.08,
        pi_hyper=0.2,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2_with_error",
        optimizer="Adamax",
        lag1=50,
        lag2=500,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1-s{nstrain}-seed{seed}.benchmark"
    resources:
        pmem=resource_calculator(data=20, nstrain=1, agg=math.prod),
        gpu_mem_mb=resource_calculator(data=20, nstrain=1, agg=math.prod),
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.history} {output.world}
        python3 -m sfacts simple_fit -m {params.model_name}  \
                --verbose --device {params.device} \
                --inpath {input.data} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --optimizer {params.optimizer} --optimizer-learning-rate {params.lr} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                -s {params.nstrain} \
                --random-seed {params.seed} \
                --outpath {output.world} --history-out {output.history}
        """

use rule fit_sfacts_strategy1 as fit_sfacts_strategy2 with:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts2-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts2-s{nstrain}-seed{seed}.history.list",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=1e-6,
        rho_hyper=1e5,
        pi_hyper=1.0,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2_with_error",
        optimizer="Adamax",
        lag1=50,
        lag2=100,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts2-s{nstrain}-seed{seed}.benchmark"

use rule fit_sfacts_strategy1 as fit_sfacts_strategy4 with:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts4-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts4-s{nstrain}-seed{seed}.history.list",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=1e-6,
        rho_hyper=1e5,
        pi_hyper=1.0,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2",
        optimizer="Adamax",
        lag1=50,
        lag2=100,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts4-s{nstrain}-seed{seed}.benchmark"


use rule fit_sfacts_strategy1 as fit_sfacts_strategy3 with:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts3-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts3-s{nstrain}-seed{seed}.history.list",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=1e-6,
        rho_hyper=1e5,
        pi_hyper=1.0,
        seed=lambda w: int(w.seed),
        model_name="simple",
        optimizer="Adamax",
        lag1=50,
        lag2=100,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts3-s{nstrain}-seed{seed}.benchmark"

use rule fit_sfacts_strategy1 as fit_sfacts_strategy5 with:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts5-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts5-s{nstrain}-seed{seed}.history.list",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=0.0005,
        rho_hyper=0.08,
        pi_hyper=0.2,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2_with_error",
        optimizer="Adamax",
        lag1=50,
        lag2=100,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts5-s{nstrain}-seed{seed}.benchmark"

use rule fit_sfacts_strategy1 as fit_sfacts_strategy6 with:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts6-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts6-s{nstrain}-seed{seed}.history.list",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=0.00005,
        rho_hyper=0.08,
        pi_hyper=0.3,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2_with_error",
        optimizer="Adamax",
        lag1=50,
        lag2=500,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts6-s{nstrain}-seed{seed}.benchmark"

use rule fit_sfacts_strategy1 as fit_sfacts_strategy7 with:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfacts7-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfacts7-s{nstrain}-seed{seed}.history.list",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=0.5,
        rho_hyper=1.0,
        pi_hyper=0.5,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2",
        optimizer="Adamax",
        lag1=50,
        lag2=500,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts7-s{nstrain}-seed{seed}.benchmark"


use rule fit_sfacts_strategy1 as fit_sfacts_strategyN with:
    output:
        world="{stem}.metagenotype-n{n}-g{g}.fit-sfactsN-pi{pi_hyper}-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype-n{n}-g{g}.fit-sfactsN-pi{pi_hyper}-s{nstrain}-seed{seed}.history.list",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=0.0005,
        rho_hyper=0.08,
        pi_hyper=lambda w: float(w.pi_hyper) / 100,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2_with_error",
        optimizer="Adamax",
        lag1=50,
        lag2=500,
        lr=5e-1,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfactsN-pi{pi_hyper}-s{nstrain}-seed{seed}.benchmark"



use rule fit_sfacts_strategy1 as fit_sfacts_strategy1_cpu with:
    output:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_cpu-s{nstrain}-seed{seed}.world.nc",
    params:
        device="cpu",
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=0.0005,
        rho_hyper=0.08,
        pi_hyper=0.2,
        learning_rate=1e-2,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2",
        optimizer="Adamax",
        lag1=40,
        lag2=200,
        lr=1e-3,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_cpu-s{nstrain}-seed{seed}.benchmark"


use rule fit_sfacts_strategy1 as fit_sfacts_strategy1_gpu with:
    output:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_gpu-s{nstrain}-seed{seed}.world.nc",
    params:
        device="cuda",
        nstrain=lambda w: int(w.nstrain),
        gamma_hyper=0.0005,
        rho_hyper=0.08,
        pi_hyper=0.2,
        learning_rate=1e-2,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2",
        optimizer="Adamax",
        lag1=40,
        lag2=200,
        lr=1e-3,
    benchmark:
        "{stem}.metagenotype-n{n}-g{g}.fit-sfacts1_gpu-s{nstrain}-seed{seed}.benchmark"



rule build_world_from_tsv:
    output:
        "{stem}.fit-{params}.to_sfacts.world.nc",
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
    shell:
        """
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        {input.script} {input.sim} {input.fit} {input.bench} {output}
        """

rule fit_sfacts_strategy8:
    output:
        world="{stem}.metagenotype{stemB}.fit-sfacts8-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        history="{stem}.metagenotype{stemB}.fit-sfacts8-s{nstrain}-g{nposition}-seed{seed}.history.list",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-5,
        rho_hyper=0.1,
        pi_hyper=0.2,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2",
        lag1=50,
        lag2=100,
        anneal_steps=2500,
        anneal_wait=500,
        lr=0.1,
        collapse=0.02,
        cull=0.001,
    resources:
        pmem=resource_calculator(data=2, nposition=0.002, nstrain=2, agg=math.prod),
        gpu_mem_mb=resource_calculator(data=2, nposition=0.002, nstrain=2, agg=math.prod),
        walltime_hr=12,
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.history} {output.world}
        python3 -m sfacts complex_fit -m {params.model_name}  \
                --verbose --device {params.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --anneal-hyperparameters gamma_hyper pi_hyper rho_hyper \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --refinement-hyperparameters gamma_hyper=1.0 rho_hyper=5.0 pi_hyper=0.1 \
                --optimizer-learning-rate {params.lr} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --num-positions {params.nposition} \
                --precision {params.precision} \
                --collapse {params.collapse} --cull {params.cull} \
                --random-seed {params.seed} \
                --history-out {output.history} \
                {input.data} {output.world}
        """

use rule fit_sfacts_strategy8 as fit_sfacts_strategy9 with:
    output:
        world="{stem}.metagenotype{stemB}.fit-sfacts9-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        history="{stem}.metagenotype{stemB}.fit-sfacts9-s{nstrain}-g{nposition}-seed{seed}.history.list",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-5,
        rho_hyper=0.1,
        pi_hyper=0.2,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2",
        lag1=50,
        lag2=500,
        anneal_steps=5000,
        anneal_wait=500,
        lr=0.1,
        collapse=0.02,
        cull=0.001,
    resources:
        pmem=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
        gpu_mem_mb=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
        walltime_hr=24,

use rule fit_sfacts_strategy8 as fit_sfacts_strategy10 with:
    output:
        world="{stem}.metagenotype{stemB}.fit-sfacts10-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        history="{stem}.metagenotype{stemB}.fit-sfacts10-s{nstrain}-g{nposition}-seed{seed}.history.list",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-3,
        rho_hyper=0.1,
        pi_hyper=0.2,
        seed=lambda w: int(w.seed),
        model_name="simple",
        lag1=50,
        lag2=500,
        anneal_steps=5000,
        anneal_wait=500,
        lr=0.1,
        collapse=0.02,
        cull=0.001,
    resources:
        pmem=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
        gpu_mem_mb=5000,  # resource_calculator(data=1, nposition=0.002, nstrain=2, agg=math.prod),
        walltime_hr=24,


rule fit_sfacts_strategy11:
    output:
        world="{stem}.metagenotype{stemB}.fit-sfacts11-s{nstrain}-seed{seed}.world.nc",
        history="{stem}.metagenotype{stemB}.fit-sfacts11-s{nstrain}-seed{seed}.history.list",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        precision=32,
        gamma_hyper=1e-4,
        rho_hyper=0.1,
        pi_hyper=0.2,
        alpha_hyper_scale=0.1,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=2500,
        anneal_wait=500,
        lr=0.05,
    resources:
        pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
        gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
        walltime_hr=12,
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.history} {output.world}
        python3 -m sfacts simple_fit -m {params.model_name}  \
                --verbose --device {params.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --anneal-hyperparameters gamma_hyper pi_hyper rho_hyper \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --optimizer-learning-rate {params.lr} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --history-out {output.history} \
                {input.data} {output.world}
        """

rule fit_sfacts_strategy12:
    output:
        world="{stem}.metagenotype{stemB}.fit-sfacts12-s{nstrain}-g{nposition}-seed{seed}.world.nc",
        history="{stem}.metagenotype{stemB}.fit-sfacts12-s{nstrain}-g{nposition}-seed{seed}.history.list",
    input:
        data="{stem}.metagenotype{stemB}.nc",
    params:
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-4,
        rho_hyper=0.1,
        pi_hyper=0.2,
        alpha_hyper_scale=0.1,
        seed=lambda w: int(w.seed),
        model_name="simple_ssdd2_with_error",
        lag1=50,
        lag2=100,
        anneal_steps=2500,
        anneal_wait=500,
        lr=0.05,
        collapse=0.02,
        cull=0.001,
    resources:
        pmem=resource_calculator(data=2, nstrain=2, agg=math.prod),
        gpu_mem_mb=resource_calculator(data=2, nstrain=2, agg=math.prod),
        walltime_hr=12,
    # conda:
    #     "conda/sfacts.yaml"
    shell:
        r"""
        export PYTHONPATH="/pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFacts"
        rm -rf {output.history} {output.world}
        python3 -m sfacts complex_fit -m {params.model_name}  \
                --verbose --device {params.device} \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --hyperparameters alpha_hyper_scale={params.alpha_hyper_scale} \
                --refinement-hyperparameters gamma_hyper=1.0 rho_hyper=5.0 pi_hyper=0.1 \
                --anneal-hyperparameters gamma_hyper pi_hyper rho_hyper \
                --anneal-steps {params.anneal_steps} --anneal-wait {params.anneal_wait} \
                --optimizer-learning-rate {params.lr} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                --num-strains {params.nstrain} \
                --num-positions {params.nposition} \
                --collapse {params.collapse} --cull {params.cull} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --history-out {output.history} \
                {input.data} {output.world}
        """
