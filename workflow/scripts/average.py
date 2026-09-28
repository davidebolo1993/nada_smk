#!/usr/bin/python3

import pandas as pd
import os
import sys

def main():

    '''
    calculate average, no need to split into multiple functions
    '''

    average_df = pd.DataFrame()

    for df in sys.argv[1:-1]:
        
        if os.stat(df).st_size != 0:

            average_df=average_df.add(pd.read_table(df, index_col=0), fill_value=0)

    average_df['cd_sum']/=average_df['N']
    average_df.rename(columns={'cd_sum': 'avg_cd'}, inplace=True)
    average_df.to_csv(sys.argv[-1], sep='\t')

if __name__ == '__main__':

    main()