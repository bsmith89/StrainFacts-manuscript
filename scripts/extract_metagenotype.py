#!/usr/bin/env python3

import sys
import sfacts as sf

if __name__ == "__main__":
    world_path = sys.argv[1]
    num_samples = int(sys.argv[2])
    num_positions = int(sys.argv[3])
    tsv_outpath = sys.argv[4]
    nc_outpath = sys.argv[5]

    world = sf.data.World.load(world_path)
    world_ss = world.isel(
        sample=slice(0, num_samples), position=slice(0, num_positions)
    )
    world_ss.metagenotypes.to_csv(tsv_outpath, sep="\t")
    world_ss.metagenotypes.dump(nc_outpath)
