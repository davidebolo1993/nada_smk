#!/usr/bin/python3

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import sys
import os

def aggregate(h5,label,apply_mean=False) -> (pd.DataFrame,np.ndarray):

    '''
    aggregation based on one or more labels - eventually use mean
    '''

    #split on label
    if label != ['None']:

        labels=['pheno_id'] + label
        donors, donors_df=map(list, zip(*h5.obs.groupby(labels)))

    else: #is None

        labels=['pheno_id']
        donors, donors_df=map(list, zip(*h5.obs.groupby(labels)))

    #initialize empty
    X_agg=np.zeros((len(donors), h5.shape[1]),dtype=np.float32)
    N_cells=np.zeros(len(donors),dtype=np.int32)

    #add additional layers
    X_agg_counts=np.zeros((len(donors), h5.shape[1]),dtype=np.float32)
    X_agg_raw=np.zeros((len(donors), h5.shape[1]),dtype=np.float32)

    #aggregation
    for i,df in enumerate(donors_df):

        h5_sub=h5[df.index]
        X_agg[i,:]=np.asarray(h5_sub.X.sum(0))
        X_agg_counts[i,:]=np.asarray(h5_sub.layers['counts'].sum(0))
        X_agg_raw[i,:]=np.asarray(h5_sub.raw.X.sum(0))
        N_cells[i]= h5_sub.obs['N_cells'].sum() if 'N_cells' in h5_sub.obs else h5_sub.shape[0]

    out_df=pd.DataFrame(donors, columns=labels)
    out_df['N_cells']=N_cells

    #use mean?
    if apply_mean:
    
        X_agg/=N_cells[:, np.newaxis] #increase dimension

    return out_df,X_agg,X_agg_counts,X_agg_raw


def main():

    if os.stat(sys.argv[1]).st_size != 0: #h5ad may be empty - which means empty tranche

        #h5in
        h5=sc.read_h5ad(sys.argv[1])

        if h5.obs.shape[0] == 0:

            #no passed
            open(sys.argv[2], 'w').close()
                    
        else:

            #change index
            #h5.obs['idx.new'] = h5.obs.index
            #h5.obs.set_index('idx.new', inplace=True)
            #h5.obs.index.rename('index', inplace=True)

            #out obj
            out_df,X_agg,X_agg_counts,X_agg_raw = aggregate(h5,sys.argv[3].split(','),apply_mean=False)

            #write output
            ad_out=ad.AnnData(X_agg, obs=out_df, var=h5.var)

            #add additional layers
            ad_out.layers['ac'] = X_agg_counts
            ad_out.layers['ar'] = X_agg_raw

            new_idx=ad_out.obs['pheno_id']
            
            if sys.argv[3] != str(None):

                    for x in sys.argv[3].split(','):

                        new_idx = [x + '__' + y.replace(' ', '_') for x,y in zip(new_idx,ad_out.obs[x])]

            ad_out.obs.index=new_idx
            ad_out.obs.index.rename('index', inplace=True)
            ad_out.write(sys.argv[2])

    else:

        #write empty
        open(sys.argv[2], 'w').close()


if __name__ == '__main__':

    main()