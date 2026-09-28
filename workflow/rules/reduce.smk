rule reduce_size:
    '''
    reduce size of pflog1pf by deleting some unnecessary layers - just for merging
    '''
    input:
        rules.applypflog1pf.output
    output:
        config['out_folder'] + 'process/pflog1pf_reduced/h5ad/{sample}/h5.h5ad'
    params:
        reduce_to=config['reduce_to']
    shell:
        '''
        python workflow/scripts/reduce.py {input} {output} {params.reduce_to}
        '''


rule reduce_size_to_X:
    '''
    reduce size of pflog1pf by deleting some unnecessary layers - just for merging
    '''
    input:
        rules.applypflog1pf.output
    output:
        config['out_folder'] + 'process/pflog1pf_reduced_X/h5ad/{sample}/h5.h5ad'
    params:
        reduce_to='X'
    shell:
        '''
        python workflow/scripts/reduce.py {input} {output} {params.reduce_to}
        '''


rule reduce_size_to_raw:
    '''
    reduce size of pflog1pf by deleting some unnecessary layers - just for merging
    '''
    input:
        rules.applypflog1pf.output
    output:
        config['out_folder'] + 'process/pflog1pf_reduced_raw/h5ad/{sample}/h5.h5ad'
    params:
        reduce_to='raw'
    shell:
        '''
        python workflow/scripts/reduce.py {input} {output} {params.reduce_to}
        '''


rule reduce_size_to_counts:
    '''
    reduce size of pflog1pf by deleting some unnecessary layers - just for merging
    '''
    input:
        rules.applypflog1pf.output
    output:
        config['out_folder'] + 'process/pflog1pf_reduced_counts/h5ad/{sample}/h5.h5ad'
    params:
        reduce_to='counts'
    shell:
        '''
        python workflow/scripts/reduce.py {input} {output} {params.reduce_to}
        '''


rule reduce_size_to_log1p_cp10k:
    '''
    reduce size of pflog1pf by deleting some unnecessary layers - just for merging
    '''
    input:
        rules.applypflog1pf.output
    output:
        config['out_folder'] + 'process/pflog1pf_reduced_log1p_cp10k/h5ad/{sample}/h5.h5ad'
    params:
        reduce_to='log1p_cp10k'
    shell:
        '''
        python workflow/scripts/reduce.py {input} {output} {params.reduce_to}
        '''