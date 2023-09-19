from setuptools import setup, find_packages

setup(
    name='gene_splicer',
    version='v2.3.9',
    packages=find_packages(),
    install_requires=[
        'gotoh @ git+https://github.com/cfe-lab/MiCall.git@v7.7.0#egg=gotoh&subdirectory=micall/alignment',
        'numpy==1.25.1',
        'python-Levenshtein==0.12.0',
        'pandas==2.0.2',
        'requests==2.31.0',
        'intactness-pipeline @ git+https://github.com/cfe-lab/HIVIntact.git@cfe-1.4',
        'pyyaml'
    ],
    package_data={
        '':
        ['genes_of_interest.yaml', 'annot.csv', 'hxb2.fasta', 'config.yaml']
    },
    entry_points={
        'console_scripts': [
            'gene_splicer_sample = gene_splicer.sample:main',
            'gene_splicer_run = gene_splicer.pipeline:main',
            'gene_splicer_study = gene_splicer.study_summary:main',
        ]
    })
