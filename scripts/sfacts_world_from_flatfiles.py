#!/usr/bin/env python3
# USAGE: <SCRIPT>.py {input.mgen} {input.gamma} {input.pi} {output}

import sfacts as sf
import sys


if __name__ == "__main__":
    metagenotypes = sf.data.Metagenotypes.load_from_tsv(sys.argv[1])
    genotypes = sf.data.Genotypes.load_from_tsv(sys.argv[2])
    communities = sf.data.Communities.load_from_tsv(sys.argv[3])

    world = sf.data.World.from_combined(
        metagenotypes,
        genotypes,
        communities,
    )
    world.data["p"] = world.data.communities @ world.data.genotypes
    world.data["mu"] = world.data.metagenotypes.sum("allele").mean("position")

    world.dump(sys.argv[4])
