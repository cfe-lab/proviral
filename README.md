# Proviral Sequence Pipeline
1. Run MiCall on sequences to produce `conseqs` and `contigs` csv files
2. Run `primer_finder.py` on the conseqs and contigs to produce a `.*filtered.csv` as well as optionally running hivseqinr to determine status of sequences