#!/usr/bin/python3

import os
import sys
import scanpy as sc
import pandas as pd
import scipy as sp
import numpy as np
import anndata

def read_h5ad(h5file) -> anndata._core.anndata.AnnData:

    '''
    scanpy read h5ad file
    '''

    return sc.read_h5ad(h5file)

def read_cd(table) -> dict:

    '''
    average cell depth to dict
    '''

    return pd.read_table(table, index_col=0)['avg_cd'].to_dict()

def splitby(h5ad,field) ->dict:

    '''
    split h5ad by field
    '''

    if field != str(None):

        return h5ad.obs.groupby(field).indices

    return {'all_cells': np.arange(h5ad.shape[0])}


def main():

    if os.stat(sys.argv[1]).st_size != 0: #h5ad may be empty

        h5=read_h5ad(sys.argv[1])
        tsv=read_cd(sys.argv[2])

        #on first pass, log transform, compute cell depth for the second round and add counts to layers
        first_pass = 'counts' not in h5.layers

        ##auxiliary array for normalization
        cd_aux=np.zeros(h5.shape[0],dtype=np.float32) #with this being set to np.float32 downstream should all be float32

        ##split-by
        idx=splitby(h5,sys.argv[3])

        for f,i in idx.items():

            #if f not in ['THP1', 'U937']: #exclude spike-ins if in object
            cd_aux[i] = tsv[f]

        ##counts per cell
        cd=np.asarray(h5.X.sum(1)).flatten()

        ##store raw counts
        if first_pass:
        
            h5.layers['counts']=h5.X.copy()

        ##normalization
        h5.X=sp.sparse.diags(cd_aux/cd).dot(h5.X)

        #log-transform and compute cell_depth
        if first_pass:

            h5.X = h5.X.log1p()
            n,cd=dict(),dict()

            for f,i in idx.items():

                #if f not in ['THP1', 'U937']: #exclude spike-ins if in object

                n[f]=h5[i].shape[0]
                cd[f]=h5[i].X.sum()
            
            pd.DataFrame({"field": n.keys(), "N": n.values(), "cd_sum": cd.values()}).to_csv(sys.argv[5], sep='\t', index=False)

        #write normalized        
        h5.write(sys.argv[4])

    else:

        open(sys.argv[4], 'w').close()

        if sys.argv[5] != 'None':

            open(sys.argv[5], 'w').close()


if __name__ == '__main__':

    main()
