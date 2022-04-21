#!/usr/bin/env python3
# {input.script} {params.args} {output.geno} {output.comm}

import pandas as pd
import sys
from tqdm import tqdm


def iter_split_mixture_s_output(lines):
    curr_strain_lines = []
    for this_line in lines:
        if this_line.startswith(">"):
            if curr_strain_lines:
                yield curr_strain_lines
                curr_strain_lines = []
        curr_strain_lines.append(this_line)
    if curr_strain_lines:
        yield curr_strain_lines
    else:
        yield ['>1.0\n']  # If the file is empty, act like there's one reference strain.


def parse_geno_and_rabund_from_dummy(strain_lines):
    assert strain_lines[0].startswith(">")
    rabund = float(strain_lines[0][1:])
    alt_positions = [int(line.strip().split(",")[1]) for line in strain_lines[1:]]
    return alt_positions, rabund


def geno_from_alt_positions(alt_positions, all_positions):
    geno = pd.Series({position: "alt" for position in alt_positions})
    return geno.reindex(all_positions, fill_value="ref")


if __name__ == "__main__":
    mgen_inpath = sys.argv[1]
    sample_args = sys.argv[2].split(":")
    geno_outpath = sys.argv[3]
    comm_outpath = sys.argv[4]

    mgen = (
        pd.read_table(
            mgen_inpath,
            dtype=str,
        )
        .astype(dict(sample=str, position=int, allele=str))
        .set_index(["sample", "position", "allele"])
        .squeeze()
        .astype(float)
        .astype(int)
    )

    geno_out = {}
    comm_out = {}
    for arg in tqdm(sample_args):
        sample_name, sample_path = arg.split("=")
        comm_out[sample_name] = {}
        with open(sample_path) as f:
            for i, strain_lines in enumerate(iter_split_mixture_s_output(f)):
                strain_name = f"{sample_name}-{i}"
                alt_positions, rabund = parse_geno_and_rabund_from_dummy(strain_lines)
                comm_out[sample_name][strain_name] = rabund
                geno = geno_from_alt_positions(
                    alt_positions, mgen.index.get_level_values("position").unique()
                )

                geno_out[strain_name] = geno
    geno_out = (
        pd.DataFrame(geno_out)
        .rename_axis(columns="strain")
        .unstack()
        .to_frame(name="allele")
        .assign(genotypes=1.0)
        .set_index("allele", append=True)
        .squeeze()
        .unstack("allele", fill_value=0)
        .reindex(columns=["alt", "ref"], fill_value=0)
        .stack()
    )
    comm_out = (
        pd.DataFrame(comm_out)
        .fillna(0)
        .rename_axis(index="strain", columns="sample")
        .unstack()
        .rename("communities")
    )

    geno_out.to_csv(geno_outpath, sep="\t", header=True)
    comm_out.to_csv(comm_outpath, sep="\t", header=True)
