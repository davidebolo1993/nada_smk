rule averagecd:
    '''
    divide sum by number of observations and get average count for category in all tranches
    '''
    input:
        expand(config['out_folder'] + 'process/celldepth/{sample}/cd.tsv', sample=sample)
    output:
        config['out_folder'] +  'process/average/cd.tsv'
    shell:
        '''
        python workflow/scripts/average.py {input} {output}
        '''

rule applylog1pf:
    '''
    apply first logpf
    '''
    input:
        avg=rules.averagecd.output,
        h5=rules.mergeh5.output.h5
    output:
        h5=config['out_folder'] + 'process/log1pf/h5ad/{sample}/h5.h5ad',
        cd=config['out_folder'] + 'process/log1pf/celldepth/{sample}/cd.tsv'
    params:
        norm=config['norm']
    shell:
        '''
        python workflow/scripts/logpf.py {input.h5} {input.avg} {params.norm} {output.h5} {output.cd}
        '''

rule averagelogcd:
    '''
    repeat averagecd but on log-transformed obs
    '''
    input:
        expand(config['out_folder'] + 'process/log1pf/celldepth/{sample}/cd.tsv', sample=sample)
    output:
        config['out_folder'] +  'process/log1pf/average/cd.tsv'
    shell:
        '''
        python workflow/scripts/average.py {input} {output}
        '''

rule applypflog1pf:
    '''
    apply second pf
    '''
    input:
        avg=rules.averagelogcd.output,
        h5=rules.applylog1pf.output.h5
    output:
        config['out_folder'] + 'process/pflog1pf/h5ad/{sample}/h5.h5ad'
    params:
        norm=config['norm']
    shell:
        '''
        python workflow/scripts/logpf.py {input.h5} {input.avg} {params.norm} {output} 'None'
        '''