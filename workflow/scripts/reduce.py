#!/usr/bin/python3

import scanpy as sc
import os
import sys

def main():

    if os.stat(sys.argv[1]).st_size != 0: #h5ad may be empty

        h5=sc.read_h5ad(sys.argv[1])

        to_keep=sys.argv[3].split(',')

        if 'X' not in to_keep:

            del(h5.X)

        if 'counts' not in to_keep:

            del(h5.layers['counts'])

        if 'log1p_cp10k' not in to_keep:

            del(h5.layers['log1p_cp10k'])

        if 'raw' not in to_keep:

            del(h5.raw)
        
        h5.write_h5ad(sys.argv[2])

    else:

        open(sys.argv[2], 'w').close()

if __name__ == '__main__':

    main()
