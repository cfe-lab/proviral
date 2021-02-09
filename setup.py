from setuptools import setup, find_packages
# from constants import PIPELINE_VERSION
setup(
    name='gene_splicer',
    version='v1.0.0',
    packages=find_packages(),
    install_requires=[
        'gotoh @ git+https://github.com/cfe-lab/MiCall.git@v7.7.0#egg=gotoh&subdirectory=micall/alignment', 'numpy==1.18.4',
        'python-Levenshtein==0.12.0', 'pandas==1.0.5', 'matplotlib==3.3.3',
        'requests==2.24.0', 'pyyaml', 'pytest'
    ],
    package_data={
        '':
        ['genes_of_interest.yaml', 'annot.csv', 'hxb2.fasta', 'config.yaml']
    },
    entry_points={
        'console_scripts': [
            'gene_splicer_run = gene_splicer.pipeline:main',
        ]
    })
