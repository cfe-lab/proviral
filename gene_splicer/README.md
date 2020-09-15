# Proviral Sequence Pipeline
1. Run MiCall on sequences to produce `conseqs` and `contigs` csv files
2. Run `primer_finder.py` on the conseqs and contigs to produce a `.*filtered.csv` as well as optionally running hivseqinr to determine status of sequences

## hxb2_custom.fasta
A custom version of hxb2 that has

# Installation
```bash
# Create a virtual environment called "proviral" (this is for Python >= 3.3)
python -m venv proviral
# Clone the repository (or download it as a zipfile)
git clone git@github.com:cfe-lab/proviral.git
# Enter the repository's top-level directory
cd proviral
# Install (include "-e" for dev mode)
pip install .
```

# Testing
* Pytest is used to run tests, make sure that you use the pytest installed in your `proviral` virtual environment (see above) to avoid `ModuleNotFound` errors, i.e.:

```
python -m pytest gene_splicer/tests
# AVOID "pytest gene_splicer/tests"
```