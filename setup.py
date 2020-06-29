from setuptools import setup, find_packages
# from constants import PIPELINE_VERSION
setup(
    name='gene_splicer',
    version='v1.0.0',
    packages=find_packages(),
    install_requires=[
        'gotoh @ git+https://github.com/cfe-lab/MiCall.git@v7.7.0#egg=gotoh&subdirectory=micall/alignment',
        'iva @ git+https://github.com/cfe-lab/iva.git@v1.1.1',
        'numpy==1.18.4',
        'python-Levenshtein==0.12.0',
        'pandas==1.0.5',
        'matplotlib==3.2.1',
        'requests==2.24.0'
    ],
    entry_points={
        'console_scripts': [
            'ihec_run = gene_splicer.master:main',
        ]
    }
)