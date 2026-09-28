from glob import glob

def get_files(x):

    if config['by'] == 'tranche':

        return glob(config['in_folder'] + x + '/Donor_Quantification/*/*.h5ad')

    else: ##this is by pool

        return glob(config['in_folder'] + '*/Donor_Quantification/' + x + '/*.h5ad')

rule geth5:
    '''
    create a list of all h5ad file names for each tranch in a freezed release
    '''
    input:
        lambda wildcards: get_files(wildcards.sample)
    output:
        config['out_folder'] + 'preprocess/fofn/all/{sample}/h5.fofn'
    shell:
        '''
        ls {input} > {output}
        '''

rule filterh5:
    '''
    filter all the h5ad files retrieved, put in include/exclude based on external table
    '''
    input:
        fofn=rules.geth5.output
    output:
        include=config['out_folder'] + 'preprocess/fofn/filtered/include/{sample}/h5.fofn',
        exclude=config['out_folder'] + 'preprocess/fofn/filtered/exclude/{sample}/h5.fofn'
    params:
        keep_file=config['keep'] if config['keep'] is not None else "None"
    shell:
        '''
        python workflow/scripts/filter.py {input.fofn} {params.keep_file} {output.include} {output.exclude}
        '''