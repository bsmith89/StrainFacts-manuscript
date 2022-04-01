use rule start_jupyter as start_jupyter_sfinder with:
    conda:
        "conda/sfinder.yaml"


use rule start_ipython as start_ipython_sfinder with:
    conda:
        "conda/sfinder.yaml"


use rule start_shell as start_shell_sfinder with:
    conda:
        "conda/sfinder.yaml"


rule metagenotype_tsv_to_sfinder_aln:
    output:
        cpickle="{stem}.metagenotype-n{n}-g{g}.sfinder.aln.cpickle",
        indexes="{stem}.metagenotype-n{n}-g{g}.sfinder.aln.indexes.txt",
    input:
        script="scripts/metagenotype_to_sfinder_alignment.py",
        data="{stem}.metagenotype-n{n}-g{g}.tsv",
    conda:
        "conda/sfinder.yaml"
    shell:
        """
        {input.script} {input.data} {output}
        """


localrules:
    metagenotype_tsv_to_sfinder_aln,


rule fit_sfinder:
    output:
        "{stem}.fit-sfinder-s{nstrain}-seed{seed}.em.cpickle",
    input:
        "{stem}.sfinder.aln.cpickle",
    conda:
        "conda/sfinder.yaml"
    params:
        nstrain=lambda w: int(w.nstrain),
        seed=lambda w: int(w.seed),
        max_runtime_s=72000,
    benchmark:
        "{stem}.fit-sfinder-s{nstrain}-seed{seed}.benchmark"
    resources:
        walltime_sec=72000,
    shell:
        """
        rm -rf {output}
        /pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFinder/StrainFinder.py \
                --force_update --merge_out --msg \
                --aln {input} \
                --seed {params.seed} \
                -N {params.nstrain} \
                --max_reps 1 --dtol 1 --ntol 2 --max_time {resources.walltime_sec} --n_keep 5 --converge \
                --em_out {output}
        """

rule fit_sfinder_exhaustive:
    output:
        "{stem}.fit-sfinder_ex-s{nstrain}-seed{seed}.em.cpickle",
    input:
        "{stem}.sfinder.aln.cpickle",
    conda:
        "conda/sfinder.yaml"
    params:
        nstrain=lambda w: int(w.nstrain),
        seed=lambda w: int(w.seed),
        max_runtime_s=72000,
    benchmark:
        "{stem}.fit-sfinder_ex-s{nstrain}-seed{seed}.benchmark"
    resources:
        walltime_sec=72000,
    shell:
        """
        rm -rf {output}
        /pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFinder/StrainFinder.py \
                --force_update --merge_out --msg \
                --aln {input} \
                --exhaustive \
                --seed {params.seed} \
                -N {params.nstrain} \
                --max_reps 1 --dtol 1 --ntol 2 --max_time {resources.walltime_sec} --n_keep 5 --converge \
                --em_out {output}
        """


rule fit_sfinder_timeit:
    output:
        time="{stem}.fit-sfinder_timeit-s{nstrain}-seed{seed}.time",
    input:
        "{stem}.sfinder.aln.cpickle",
    conda:
        "conda/sfinder.yaml"
    params:
        nstrain=lambda w: int(w.nstrain),
        seed=lambda w: int(w.seed),
        out="{stem}.fit-sfinder_timeit-s{nstrain}-seed{seed}.em.cpickle",
    resources:
        walltime_sec=172_800,
    benchmark:
        "{stem}.fit-sfinder_timeit-s{nstrain}-seed{seed}.benchmark"
    shell:
        """
        rm -rf {output}
        `which time` -o {output} -- /pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFinder/StrainFinder.py \
                --force_update --merge_out --msg \
                --aln {input} \
                --seed {params.seed} \
                --random \
                -N {params.nstrain} \
                --max_reps 1 --max_time 1 --n_keep 1 \
                --em_out {params.out}
        """


rule fit_sfinder_global:
    output:
        "{stem}.fit-sfinder2-s{nstrain}-seed{seed}.em.cpickle",
    input:
        "{stem}.sfinder.aln.cpickle",
    conda:
        "conda/sfinder.yaml"
    params:
        nstrain=lambda w: int(w.nstrain),
        seed=lambda w: int(w.seed),
    benchmark:
        "{stem}.fit-sfinder2-s{nstrain}-seed{seed}.benchmark"
    shell:
        """
        rm -rf {output}
        /pollard/home/bsmith/Projects/haplo-benchmark/include/StrainFinder/StrainFinder.py \
                --force_update --merge_out --msg \
                --aln {input} \
                --seed {params.seed} \
                -N {params.nstrain} \
                --max_reps 10 --dtol 1 --ntol 2 --max_time 7200 --n_keep 3 --converge \
                --em_out {output}
        """


rule parse_sfinder_cpickle:
    output:
        pi="{stem}.fit-sfinder{params}.pi.tsv",
        gamma="{stem}.fit-sfinder{params}.gamma.tsv",
    input:
        script="scripts/strainfinder_result_to_flatfiles.py",
        cpickle="{stem}.fit-sfinder{params}.em.cpickle",
        indexes="{stem}.sfinder.aln.indexes.txt",
    conda:
        "conda/sfinder.yaml"
    shell:
        """
        {input.script} {input.cpickle} {input.indexes} {output.pi} {output.gamma}
        """


localrules:
    parse_sfinder_cpickle,
