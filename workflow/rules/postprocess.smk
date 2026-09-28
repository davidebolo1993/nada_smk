rule aggregate:
    input:
        rules.applypflog1pf.output
    output:
        config['out_folder'] + 'postprocess/aggregate/h5ad/{sample}/h5.h5ad',
    params:
        aggregation=config['aggregate']
    shell:
        '''
        python workflow/scripts/aggregate.py {input} {output} {params.aggregation}
        '''

rule final_aggregate:
    input:
        expand(config['out_folder'] + 'postprocess/aggregate/h5ad/{sample}/h5.h5ad', sample=sample)
    output:
        split=directory(config['out_folder'] + 'postprocess/aggregate/results/split'),
        h5avg=config['out_folder'] + 'postprocess/aggregate/results/cardinal_all_tranches.agg_mean.h5ad',
        tsvavg=config['out_folder'] + 'postprocess/aggregate/results/cardinal_all_tranches.agg_mean.tsv',
        statsavg=config['out_folder'] + 'postprocess/aggregate/results/cardinal_all_tranches.agg_mean.stats.tsv',
        covaravg=config['out_folder'] + 'postprocess/aggregate/results/cardinal_all_tranches.agg_mean.covar.tsv'
    threads:
        3
    params:
        aggregation=config['aggregate']
    shell:
        '''
        python workflow/scripts/final_aggregate.py {params.aggregation} {output.h5avg} {output.tsvavg} {output.statsavg} {output.covaravg} {input} {output.split}
        '''