rule mergeh5:
    '''
    create one h5ad file for each tranche - can be empty or not but we keep track of what we have done. Parallelize also reading within tranches. Compute cell-depth per category
    '''
    input:
        rules.filterh5.output.include
    output:
        h5=config['out_folder'] + 'merge/h5ad/{sample}/h5.h5ad',
        cd=config['out_folder'] + 'process/celldepth/{sample}/cd.tsv'
    params:
        fields=config['add'],
        norm=config['norm'],
        retain=config['filter'],
        include=config['include_info']
    threads:
        8
    shell:
        '''
        python workflow/scripts/mergeh5.py {input} {output.h5} "{params.fields}" {threads} {output.cd} {params.norm} {params.retain} {params.include}
        '''
