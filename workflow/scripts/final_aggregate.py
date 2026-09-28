#!/usr/bin/python3

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.decomposition import PCA
import sys
import os
import multiprocessing


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
        #left untouched
        X_agg[i,:]=np.asarray(h5_sub.X.sum(0))
        #addition for sum counts
        X_agg_counts[i,:]=np.asarray(h5_sub.layers['ac'].sum(0))
        X_agg_raw[i,:]=np.asarray(h5_sub.layers['ar'].sum(0))
        N_cells[i]= h5_sub.obs['N_cells'].sum() if 'N_cells' in h5_sub.obs else h5_sub.shape[0]

    out_df=pd.DataFrame(donors, columns=labels)
    out_df['N_cells']=N_cells

    #use mean?
    if apply_mean:
    
        X_agg/=N_cells[:, np.newaxis] #increase dimension

    return out_df,X_agg,X_agg_counts,X_agg_raw


def save_stats(h5,file):

    '''
    save statistics
    '''

    genes=h5.var.index.tolist()
    mean=h5.X.mean(0)
    var=h5.X.var(0)

    pd.DataFrame({'gene': genes, 'mean': mean, 'var': var}).to_csv(file, sep='\t', index=False)


def save_covar(h5,components,file):

    '''
    save covariates
    '''

    if components>0:

        N=min(components, *h5.shape)
        pd.DataFrame(PCA(N).fit_transform(h5.X), index=h5.obs.index, columns=[f'pca_{i}' for i in range(N)]).to_csv(file, sep='\t', index=False)

def chunks(l,n) -> list:

    '''
    split list in chunks based on number of threads
    '''

    return [l[i:i+n] for i in range(0, len(l), n)]


def save_donor_ncells(h5, file):

    '''
    save number of cells per donor to a dedicated tsv file    
    '''
    
    donor=h5.obs['pheno_id']
    N_cells=h5.obs['N_cells']
    df=pd.DataFrame([x.split('__') for x in donor], columns=['tranche.id', 'pool.id', 'donor.id'])
    df['unique.id'] = [str(x) for x in donor]
    df['N_cells'] = [str(x) for x in N_cells]
    df = df[['unique.id', 'tranche.id', 'pool.id', 'donor.id', 'N_cells']]
    df.to_csv(file, sep='\t', index=False)


def save_gene_counts(h5,file):

    '''
    save gene counts
    '''

    genes=h5.var.index.tolist()
    tot_counts=h5.layers['ac'].sum(0)
    tot_counts_raw=h5.layers['ar'].sum(0)

    pd.DataFrame({'gene': genes, 'sum_counts': tot_counts, 'sum_counts_raw': tot_counts_raw}).to_csv(file, sep='\t', index=False)


def process(h5,sli):

    '''
    multiprocessing and writing
    '''

    label=sli[0]
    out_df,X_agg,X_agg_counts,X_agg_raw= aggregate(h5, [label], apply_mean=True) #mean is only applied to X_agg 

    ad_tmp=ad.AnnData(X_agg, obs=out_df, var=h5.var)
    ad_tmp.obs.index=ad_tmp.obs['pheno_id'] + '__' + ad_tmp.obs[label]
    ad_tmp.obs.index.rename('index', inplace=True)
    ad_tmp.layers['ac'] = X_agg_counts
    ad_tmp.layers['ar'] = X_agg_raw

    #create output split folder
    os.makedirs(sys.argv[-1] + '/' + label, exist_ok=True)

    #splits
    splits, splits_df=map(list, zip(*ad_tmp.obs.groupby(label)))

    for s,d in zip(splits,splits_df):

        s=s.replace(' ', '_')
        s=s.replace('/', '_')
        adata_sub=ad_tmp[d.index].copy()
        adata_sub.obs.index=adata_sub.obs['pheno_id']
        adata_sub.obs.index.rename('index', inplace=True)

        save_donor_ncells(adata_sub,sys.argv[-1] + '/' + label + '/' + s + '.ncells_x_donor.tsv')
        save_gene_counts(adata_sub,sys.argv[-1] + '/' + label + '/' + s + '.counts_x_gene.tsv')

        adata_sub.write(sys.argv[-1] + '/' + label + '/' + s + '.agg_mean.h5ad')
        adata_sub.to_df().T.to_csv(sys.argv[-1] + '/' + label + '/' + s + '.agg_mean.tsv', sep='\t')
        #output sum - ac and ar where not touched by average, they are just sums
        adata_sub.to_df('ac').T.to_csv(sys.argv[-1] + '/' + label + '/' + s + '.count.agg_sum.tsv', sep='\t')
        adata_sub.to_df('ar').T.to_csv(sys.argv[-1] + '/' + label + '/' + s + '.raw.agg_sum.tsv', sep='\t')

        save_stats(adata_sub,sys.argv[-1] + '/' + label + '/' + s + '.agg_mean.stats.csv')
        save_covar(adata_sub,10,sys.argv[-1] + '/' + label + '/' + s + '.agg_mean.covar.tsv')

def main():

    #read input
    h5l=list()

    for h5 in sys.argv[6:-1]:
        
        if os.stat(h5).st_size != 0:

            h5l.append(sc.read_h5ad(h5))

    #concatenate
    h5=ad.concat(h5l)

    #no need to check if previous is empty or not
    #calculate average for a single donor using all the N_cells for that donor
    out_df,X_agg,_,_ = aggregate(h5,['None'],apply_mean=True)

    #merged object
    ad_out=ad.AnnData(X_agg, obs=out_df, var=h5.var)
    ad_out.obs.index=ad_out.obs['pheno_id']
    ad_out.obs.index.rename('index', inplace=True)
    
    #save h5ad
    ad_out.write(sys.argv[2])
    #save counts tsv
    ad_out.to_df().T.to_csv(sys.argv[3], sep='\t')
    #stats to file
    save_stats(ad_out,sys.argv[4])
    #covar to file
    save_covar(ad_out,10,sys.argv[5])

    #create output folder
    os.makedirs(sys.argv[-1], exist_ok=True)

    #if None nothing to do, otherwise split
    if sys.argv[1] != str(None):

        labels=sys.argv[1].split(',')
        chunk_size=1
        slices=chunks(labels,chunk_size)
        processes = []

        for sli in slices:
           
            p=multiprocessing.Process(target=process, args=(h5,sli))    
            p.start()
            processes.append(p)

        for p in processes:

            p.join()


if __name__ == '__main__':

    main()