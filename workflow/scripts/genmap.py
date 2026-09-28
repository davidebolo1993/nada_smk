#!/usr/bin/python3

import pandas as pd
import sys

def read_fofn(fofn) -> list:

    '''
    list of available reports
    '''

    files=list()

    with open(fofn, 'r') as fin:

        for line in fin:

            files.append(line.rstrip())

    return files


def read_report(table)-> pd.DataFrame:

    '''
    read report and force to string certain fields
    '''

    return pd.read_table(table, dtype={"Pool_ID.Donor_Id": str, "Experiment ID": str, "Pool ID": str, "Donor id": str})


def read_sample(table) -> pd.DataFrame:

    '''
    read sample file
    '''

    df=pd.read_csv(table, sep=' ', dtype={"ID_2": str, "sex": str}).rename(columns={"ID_2":"geno_id"}).iloc[1:]
    df['geno_id_short'] = df['geno_id'].apply(lambda x: x.split("_")[0])
    return df


def concat_table(table_list) -> pd.DataFrame:

    '''
    concatenate list of tables
    '''

    return pd.concat(table_list)


def main():

    fofn=read_fofn(sys.argv[1])

    reports=[]

    for f in fofn:

        reports.append(read_report(f))

    report=concat_table(reports)
    report['pheno_id']=report['Experiment ID'] + '__' + report['Pool ID'] + '__' + report['Pool_ID.Donor_Id'].apply(lambda x: x.split('_')[-1])
    ukbb_samples=read_sample(sys.argv[2])
    elgh_samples=read_sample(sys.argv[3])
    all_samples=concat_table([ukbb_samples,elgh_samples])
    
    out=pd.merge(report, all_samples, how='inner', left_on='Donor id', right_on='geno_id_short')
    out[['geno_id', 'pheno_id', 'cohort']].to_csv(sys.argv[4], header=False, index=False, sep='\t')
    out.loc[out.cohort=='GT_UKBB'][['geno_id', 'pheno_id']].to_csv(sys.argv[5], header=False, index=False, sep='\t')
    out.loc[out.cohort=='GT_ELGH'][['geno_id', 'pheno_id']].to_csv(sys.argv[6], header=False, index=False, sep='\t')


if __name__ == '__main__':

    main()
