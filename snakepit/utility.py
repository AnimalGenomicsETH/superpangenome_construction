def get_reference_ID():
    for sample, preset in config['pangenome_samples'].items():
        if preset == 'reference':
            return sample
