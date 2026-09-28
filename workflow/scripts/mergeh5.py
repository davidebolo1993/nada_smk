#!/usr/bin/python3

import pandas as pd
import scanpy as sc
import anndata
import os
import sys
import multiprocessing
import math

def read_fofn(fofn) -> list:

    '''
    tranche-level or pool-level h5ad fofn
    '''

    files=list()

    with open(fofn, 'r') as fin:

        for line in fin:

            files.append(line.rstrip())

    return files


def read_fields(fstring) -> list:

    '''
    fields to use
    '''
    
    if fstring == str(None):

        return list()

    return fstring.split(',')


def read_h5ad(h5file) -> anndata._core.anndata.AnnData:

    '''
    scapy read h5ad file
    '''

    return sc.read_h5ad(h5file)


def read_table(table, fields) -> pd.DataFrame:

    '''
    read report table if necessary
    '''

    if os.path.exists(table):

        if len(fields) != 0:

            if 'Pool_ID.Donor_Id' in fields and len(fields) > 1:

                return pd.read_table(table, usecols=fields).set_index('Pool_ID.Donor_Id', drop=False)
            
            else:

                raise ValueError('if add is not None, it must contain at least Pool_ID.Donor_Id and another valid column')

    return pd.DataFrame()


def read_filter(table) -> pd.DataFrame:

    '''
    read filter
    '''

    if table == str(None):

        return pd.DataFrame()

    return pd.read_table(table,header=None).set_index(0, drop=False)


def concat_table(table_list) -> pd.DataFrame:

    '''
    concatenate list of tables
    '''

    if len(table_list) == 0:

        return pd.DataFrame()
    
    return pd.concat(table_list)


def get_info(file,fields) -> pd.DataFrame:

    '''
    get tables having information for current file
    '''

    tables=[]
    tranch=os.path.dirname(os.path.dirname(os.path.dirname(file)))
    tranch_id=os.path.basename(tranch)

    #ukbb table
    ukbb_path=os.path.join(tranch,'Summary_plots/Summary/UKBB_REPORT/' + tranch_id + '_UKBB_Report.tsv')

    if os.path.exists(ukbb_path):

        ukbb_table=read_table(ukbb_path,fields)
        tables.append(ukbb_table)

    #elgh_table
    elgh_path=os.path.join(tranch,'Summary_plots/Summary/ELGH_REPORT/' + tranch_id + '_ELGH_Report.tsv')

    if os.path.exists(elgh_path):

        elgh_table=read_table(elgh_path,fields)
        tables.append(elgh_table)

    return concat_table(tables)


def get_id(filename) -> str:

    '''
    get id
    '''

    fname=os.path.basename(filename)
    ffields=fname.split('.')

    return ffields[0] + '_' + ffields[1]


def chunks(l,n) -> list:

    '''
    split list in chunks based on number of threads
    '''

    return [l[i:i+n] for i in range(0, len(l), n)]


def process(sli,info,filter_table,allh5,keep_cols):

    '''
    multiprocessing reading and adding info
    '''

    for s in sli:

        h5=sc.read_h5ad(s)

        id_=get_id(s)

        if info.shape[0] != 0:

            id_info=info[info.index==id_].head(1)

            for inf in info.columns:

                var=str(id_info[inf].tolist()[0])

                if inf == 'Vacutainer ID' and var.startswith('000'):

                    inf=inf.lstrip('0')
                        
                h5.obs[inf] = var

        else:

            h5.obs['Pool_ID.Donor_Id'] = id_


        if keep_cols != ['all']:

            #keep some
            keep_cols=list(set(keep_cols + ['tranche.id', 'experiment.id','Pool_ID.Donor_Id'])) #those 3 are required
            h5.obs=h5.obs[keep_cols]

        #build unique identifier for downstream eqtl analysis
        h5.obs['pheno_id'] = h5.obs['tranche.id'].astype(str) + '__' + h5.obs['experiment.id'].astype(str) + '__' + h5.obs['Pool_ID.Donor_Id'].astype(str).apply(lambda x: x.split('_')[-1]) 
        
        #new index
        pool_info=[h.rsplit('_', 1)[0] for h in h5.obs['Pool_ID.Donor_Id']]
        h5.obs.index=[x + '--' + y + '--' + z for x,y,z in zip(h5.obs['tranche.id'],pool_info,h5.obs.index)]

        if filter_table.shape[0] == 0:

            h5.obs['Is_passed'] = True
            allh5.append(h5)
        
        else:

            passes=[True if x in filter_table.index else False for x in h5.obs.index]
            h5.obs['Is_passed']=passes
            newh5=h5[h5.obs['Is_passed']].copy()

            if filter_table.shape[1] == 2: #then we also have aggregation information to add

                aggregation_label=[filter_table.loc[x,1] for x in newh5.obs.index]
                newh5.obs['aggregation_label'] = aggregation_label

            allh5.append(newh5)


def compute_cell_depth(h5ad,field) -> (dict,dict):

    '''
    Compute cell depth
    '''
    n,cd=dict(),dict() #number of cells, cell_depth - by field

    if field != str(None):

        idx=h5ad.obs.groupby(field).indices

        for f,i in idx.items():

            #if f not in ['THP1', 'U937']: #exclude spike-ins if in object

            n[f]=h5ad[i].shape[0]
            cd[f]=h5ad[i].X.sum()

    else:

        f='all_cells'
        n[f] = h5ad.shape[0]
        cd[f] = h5ad.X.sum()

    return n,cd

def main():

    files=read_fofn(sys.argv[1])
    
    # all the files from the same tranch have same table, we can just read it once
    if len(files) != 0 :

        info=get_info(files[0],read_fields(sys.argv[3]))
        filter_table=read_filter(sys.argv[7])       
        manager = multiprocessing.Manager()
        allh5=manager.list()
        chunk_size=len(files)/int(sys.argv[4])
        slices=chunks(files,math.ceil(chunk_size))
        processes = []

        for sli in slices:

            p=multiprocessing.Process(target=process, args=(sli,info,filter_table,allh5, sys.argv[8].split(',')))
            p.start()
            processes.append(p)

        for p in processes:

            p.join()

        #concat
        finalh5=anndata.concat(allh5)

        #write
        finalh5.write(sys.argv[2])

        #compute cell_depth
        n,cd=compute_cell_depth(finalh5,sys.argv[6])
        #store to file
        pd.DataFrame({"field": n.keys(), "N": n.values(), "cd_sum": cd.values()}).to_csv(sys.argv[5], sep='\t', index=False)

    else:

        open(sys.argv[2], 'w').close()
        open(sys.argv[5], 'w').close()



if __name__ == '__main__':

    main()