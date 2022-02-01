#!/usr/bin/env python3

import sys
from sfacts import Metagenotypes

if __name__ == "__main__":
    sizes = Metagenotypes.peak_netcdf_sizes(sys.argv[1])
    for k in sizes:
        print(k, sizes[k], sep='\t')
