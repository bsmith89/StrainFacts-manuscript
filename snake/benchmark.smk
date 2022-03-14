rule run_sfinder_fit_benchmark_matrix:
    output:
        touch("data/benchmark_fit_matrix.sfinder.flag"),
    input:
        fit=[
            (
                "data/sfacts_simulate-model_simplest_simulation-n{n}-g250-s{sim_s}"
                "-pi40-mu100-eps10-seed{sim_seed}"
                ".metagenotype-n{n}-g250"
                ".fit-sfinder-s{fit_s}-seed{fit_seed}"
                ".{suffix}"
            ).format(
                sim_s=sim_s,
                sim_seed=sim_seed,
                fit_seed=fit_seed,
                n=sim_s * 5,
                fit_s=int(sim_s * fit_s_ratio),
                suffix=suffix,
            )
            for sim_s, fit_s_ratio, sim_seed, fit_seed, suffix in product(
                [10, 20, 40, 80],
                [1.0, 1.5],
                range(5),
                range(5),
                ['evaluation.tsv', 'benchmark'],
            )
        ],


rule run_sfacts_cpu_fit_benchmark_matrix:
    output:
        touch("data/benchmark_fit_matrix.sfacts_cpu.flag"),
    input:
        fit=[
            (
                "data/sfacts_simulate-model_simplest_simulation-n{n}-g250-s{sim_s}"
                "-pi40-mu100-eps10-seed{sim_seed}"
                ".metagenotype-n{n}-g250"
                ".fit-sfacts44_cpu-s{fit_s}-seed{fit_seed}"
                ".{suffix}"
            ).format(
                sim_s=sim_s,
                sim_seed=sim_seed,
                fit_seed=fit_seed,
                n=sim_s * 5,
                fit_s=int(sim_s * fit_s_ratio),
                suffix=suffix,
            )
            for sim_s, fit_s_ratio, sim_seed, fit_seed, suffix in product(
                [10, 20, 40, 80, 200],
                [1.0, 1.5],
                range(5),
                range(5),
                ['evaluation.tsv', 'benchmark'],
            )
        ],


rule run_sfacts_gpu_fit_benchmark_matrix:
    output:
        touch("data/benchmark_fit_matrix.sfacts_gpu.flag"),
    input:
        fit=[
            (
                "data/sfacts_simulate-model_simplest_simulation-n{n}-g250-s{sim_s}"
                "-pi40-mu100-eps10-seed{sim_seed}"
                ".metagenotype-n{n}-g250"
                ".fit-sfacts44_gpu-s{fit_s}-seed{fit_seed}"
                ".{suffix}"
            ).format(
                sim_s=sim_s,
                sim_seed=sim_seed,
                fit_seed=fit_seed,
                n=sim_s * 5,
                fit_s=int(sim_s * fit_s_ratio),
                suffix=suffix,
            )
            for sim_s, fit_s_ratio, sim_seed, fit_seed, suffix in product(
                [10, 20, 40, 80, 200],
                [1.0, 1.5],
                range(5),
                range(5),
                ['evaluation.tsv', 'benchmark'],
            )
        ],


# rule run_sfacts_big_fit_benchmark_matrix:
#     output:
#         touch("data/benchmark_fit_matrix.sfacts_big.flag"),
#     input:
#         fit=[
#             (
#                 "data/sfacts_simulate-model_simplest_simulation-n{n}-g{g}-s{sim_s}-pi40-mu100-eps10-seed{sim_seed}"
#                 ".metagenotype-n{n}-g{g}"
#                 ".fit-{fit_type}-s{fit_s}-seed{fit_seed}"
#                 ".evaluation.tsv"
#             ).format(
#                 fit_type=fit_type,
#                 sim_s=sim_s,
#                 sim_seed=sim_seed,
#                 fit_seed=fit_seed,
#                 n=sim_s * 5,
#                 g=g,
#                 fit_s=int(sim_s * fit_s_ratio),
#             )
#             for fit_type, sim_s, fit_s_ratio, g, sim_seed, fit_seed in product(
#                 ['sfacts41_gpu', 'sfacts41_big', 'sfacts45_big', 'sfacts44_big', 'sfacts46_big'],
#                 [40, 80, 200, 500],
#                 [1.0, 1.5],
#                 [250, 1000],
#                 range(2),
#                 range(2),
#             )
#         ],


rule run_sfinder_mem_benchmark_matrix:
    output:
        touch("data/benchmark_mem_matrix.sfinder.flag"),
    input:
        mem=[
            (
                "data/sfacts_simulate-model_simplest_simulation-n{n}-g{g}-s{s}"
                "-pi40-mu100-eps10-seed{sim_seed}"
                ".metagenotype-n{n}-g{g}"
                ".fit-sfinder_timeit-s{s}-seed{fit_seed}"
                ".time"
            ).format(g=g, s=s, n=n, fit_seed=fit_seed, sim_seed=sim_seed)
            for n, g, s, fit_seed, sim_seed in product(
                [100, 200, 500],
                [250, 500, 1000],
                [20, 40, 100, 200],
                range(3),
                range(3),
            )
        ],


rule run_sfacts_cpu_mem_benchmark_matrix:
    output:
        touch("data/benchmark_mem_matrix.sfacts_cpu.flag"),
    input:
        mem=[
            (
                "data/sfacts_simulate-model_simplest_simulation-n{n}-g{g}-s{s}"
                "-pi40-mu100-eps10-seed{sim_seed}"
                ".metagenotype-n{n}-g{g}"
                ".fit-sfacts44_timeit-s{s}-seed{fit_seed}"
                ".time"
            ).format(g=g, s=s, n=n, fit_seed=fit_seed, sim_seed=sim_seed)
            for n, g, s, fit_seed, sim_seed in product(
                [100, 200, 500, 1000, 2500, 10000],
                [250, 500, 1000],
                [20, 40, 100, 200, 400],
                range(3),
                range(3),
            )
        ],


# FIXME: Do larger N and fix the n_to_s_ratio
rule run_sfacts_gpu_mem_benchmark_matrix:
    output:
        touch("data/benchmark_mem_matrix.sfacts_gpu.flag"),
    input:
        memory=[
            (
                "data/sfacts_simulate-model_simplest_simulation-n{n}-g{g}-s{s}"
                "-pi40-mu100-eps10-seed{sim_seed}"
                ".metagenotype-n{n}-g{g}"
                ".fit-sfacts44_gpumem-s{s}-seed{fit_seed}"
                ".gpumem"
            ).format(g=g, s=s, n=n, fit_seed=fit_seed, sim_seed=sim_seed)
            for g, s, n, fit_seed, sim_seed in product(
                [250, 500, 1000],
                [20, 40, 100, 200],
                [100, 200, 500],
                range(3),
                range(3),
            )
        ],


rule all_cpu_benchmarks:
    output:
        touch("data/benchmark_all_cpu_matrix.flag"),
    input:
        sfinder_fit="data/benchmark_fit_matrix.sfinder.flag",
        sfacts_fit="data/benchmark_fit_matrix.sfacts_cpu.flag",
        sfinder_mem="data/benchmark_mem_matrix.sfinder.flag",
        sfacts_mem="data/benchmark_mem_matrix.sfacts_cpu.flag",


rule all_gpu_benchmarks:
    output:
        touch("data/benchmark_all_gpu_matrix.flag"),
    input:
        fit="data/benchmark_fit_matrix.sfacts_gpu.flag",
        mem="data/benchmark_mem_matrix.sfacts_gpu.flag",


rule all_benchmarks:
    output:
        touch("data/benchmark_all_matrix.flag"),
    input:
        cpu="data/benchmark_all_cpu_matrix.flag",
        gpu="data/benchmark_all_gpu_matrix.flag",
