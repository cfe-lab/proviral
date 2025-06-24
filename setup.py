from setuptools import setup, find_packages

setup(
    name='cfeproviral',
    version='v2.4.4',
    classifiers=[
        'License :: OSI Approved :: MIT License',
    ],
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'gotoh @ git+https://github.com/cfe-lab/MiCall.git@v7.7.0#egg=gotoh&subdirectory=micall/alignment',
        'numpy==1.25.1',
        'python-Levenshtein==0.12.0',
        'pandas==2.2.2',
        'requests==2.32.4',
        'cfeintact @ git+https://github.com/cfe-lab/CFEIntact.git@v1.23.2',
        'pyyaml'
    ],
    package_data={
        '':
        ['genes_of_interest.yaml', 'annot.csv', 'hxb2.fasta', 'config.yaml']
    },
    entry_points={
        'console_scripts': [
            'cfeproviral = cfeproviral.__main__:entry',
        ]
    })
