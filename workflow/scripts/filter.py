#!/usr/bin/python3

import sys
import re

def read_fofn(fofn) -> list:

    '''
    tranche-level or pool-level h5ad fofn
    '''

    files=list()

    with open(fofn, 'r') as fin:

        for line in fin:

            files.append(line.rstrip())

    return files


def read_filter(filter_file) -> list:

    '''
    subset to specific files
    '''

    files=list()

    if filter_file == str(None):

        return files

    #file can be anything, a table or something different - maybe another fofn - look for valid paths

    with open(filter_file, 'r') as ff:

        for line in ff:

           m = re.search('/(.+?).h5ad', line.rstrip())

           if m:

            files.append('/' + m.group(1) + '.h5ad')

    return files 


def write_fofn(outlist, outfile):

    '''
    write out fofn
    '''

    ##create empty file
    if len(outlist) == 0:

        open(outfile, 'w').close()
    
    else:

        with open(outfile, 'w') as fout:

            for line in outlist:

                fout.write(line + '\n')


def get_common(fofnlist,fflist) -> list:

    '''
    get common elements
    '''

    return list(set(fofnlist) & set(fflist))


def get_diff(fofnlist,fflist) -> list:

    '''
    get difference
    '''

    return list(set(fofnlist) - set(fflist))

    
def main():

    fofn=read_fofn(sys.argv[1])
    ff=read_filter(sys.argv[2])

    if len(ff) == 0:

        write_fofn(fofn,sys.argv[3])
        write_fofn(list(),sys.argv[4])

    else:

        intersection=get_common(fofn,ff)
        write_fofn(intersection,sys.argv[3])
        difference=get_diff(fofn,ff)
        write_fofn(difference,sys.argv[4])

                
if __name__ == '__main__':

    main()