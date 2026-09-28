rule prepare_mapping:
    input:
        fofn=config['reports'],
        ukbb_tsv=config['ukbb_samples'],
        elgh_tsv=config['elgh_samples']
    output:
        whole='resources/mapping/all.tsv',
        ukbb='resources/mapping/ukbb.tsv',
        elgh='resources/mapping/elgh.tsv'
    shell:
        '''
        python workflow/scripts/genmap.py {input.fofn} {input.ukbb_tsv} {input.elgh_tsv} {output.whole} {output.ukbb} {output.elgh}
        '''