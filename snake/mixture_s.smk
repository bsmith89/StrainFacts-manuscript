rule start_shell_mixture_s:
    container:
        "docker://bsmith89/mixture_s:d60eada44d4a5222a4729565071d084773abf066"
    shell:
        """
        export PYTHONPATH="$PWD/include/MixtureS"
        bash
        """


rule construct_dummy_fasta:
    output:
        "data/dummy_reference_fasta.g{g}.fn",
    params:
        name="dummy",
        g=lambda w: int(w.g),
        per_line=80,
        base="A",
    run:
        with open(output[0], "w") as f:
            print(f">{params.name}", file=f)
            for line_i in range(params.g // params.per_line):
                print(f"{params.base}" * params.per_line, file=f)
            print(f"{params.base}" * (params.g % params.per_line), file=f)


rule metagenotype_to_mixture_s_input:
    output:
        "{stem}.metagenotype-n{n}-g{g}.sample-{sample}.mixtureS_input.tsv",
    input:
        data="{stem}.metagenotype-n{n}-g{g}.tsv",
    params:
        sample=lambda w: w.sample,
    run:
        data = (
            pd.read_table(
                input[0],
                dtype=str,
            )
            .astype(dict(sample=str, position=int, allele=str))
            .set_index(["sample", "position", "allele"])
            .squeeze()
            .astype(float)
            .astype(int)
            .xs(params.sample, level="sample")
            .unstack("allele")
            .rename(columns={"ref": "A", "alt": "T"})
            .assign(G=0, C=0)
            .sort_index(axis="columns")
            .sort_index(axis="index")
        )
        # FIXME: No-op for debugging:
        poly = data.apply(lambda x: x / x.sum(), axis=1).max(1) < 1.1
        (data[poly].to_csv(output[0], sep="\t", header=False))


rule fit_mixture_s:
    output:
        dir=directory("{stem}.metagenotype-{params}.sample-{sample}.fit-mixtureS.d"),
        haplo="{stem}.metagenotype-{params}.sample-{sample}.fit-mixtureS.output",
    input:
        data="{stem}.metagenotype-{params}.sample-{sample}.mixtureS_input.tsv",
        genome="data/dummy_reference_fasta.g20000.fn",
    container:
        "docker://bsmith89/mixture_s:d60eada44d4a5222a4729565071d084773abf066"
    shell:
        """
        export PYTHONPATH="$PWD/include/MixtureS"
        mkdir -p {output.dir}/dummy/
        ln -rs {input.data} {output.dir}/dummy/filter_polymorphic_sites
        python3 -m mixture_model \
                --restart \
                --sample_name dummy \
                --genome_len 20000 \
                --genome_name dummy_genome \
                --genome_file_loc {input.genome} \
                --bam_file DOES_NOT_EXIST.bam \
                --res_dir {output.dir}
        cp {output.dir}/dummy/dummy_haplotypes {output.haplo}
        """


rule recombine_mixture_s_with_dummy_reference_and_fixed_number_samples:
    output:
        geno="{stem}.metagenotype-n{n}-g{g}.fit-mixtureS.gamma.tsv",
        comm="{stem}.metagenotype-n{n}-g{g}.fit-mixtureS.pi.tsv",
    input:
        script="scripts/recombined_mixtureS_outputs_from_dummy.py",
        samples=lambda w: [
            f"{w.stem}.metagenotype-n{w.n}-g{w.g}.sample-{sample}.fit-mixtureS.output"
            for sample in range(int(w.n))
        ],
        mgen="{stem}.metagenotype-n{n}-g{g}.tsv",
    params:
        args=lambda w: ":".join(
            [
                f"{sample}={w.stem}.metagenotype-n{w.n}-g{w.g}.sample-{sample}.fit-mixtureS.output"
                for sample in range(int(w.n))
            ]
        ),
    shell:
        """
        {input.script} {input.mgen} {params.args} {output.geno} {output.comm}
        """


rule construct_dummy_mixture_s_benchmark:
    output:
        "{stem}.fit-mixtureS.benchmark",
    shell:
        dd(
            """
        cat > {output} <<EOF
        s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time
        0.0	0:00:00	00.00	00.00	00.00	00.00	0.00	0.00	0.00	0.00
        EOF
        """
        )
