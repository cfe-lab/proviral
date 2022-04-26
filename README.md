## Readme

### Dependencies
1. minimap2 (https://github.com/lh3/minimap2) (must be available via commandline)
2. blast tools (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
3. R and RSCRIPT (https://www.r-project.org/)

### Singularity builds
* Build all singularity images inside of the `simages` folder

### Filtering
* At the core of the proviral pipeline, data is read from `contigs.csv` and `conseqs.csv` files produced by MiCall
* First the pipeline reads through all of the contigs, then the conseqs
* When it does this (see the `find_primers()` function) it applies the following logic in this order for filtering/tagging:
  1. If a sample is not proviral, skip it. Do not attempt to find primers or anything, just log a message saying `sample X was skipped because it was non-proviral`
  2. If a sample has 0 in the remap column of the `cascade.csv` file, tag that sequence with an error: `No contig/conseq constructed`, do not analyze it or try to find primers, and write it to the `*primer_analysis.csv` file (which records all failures)
  3. If the `consensus-percent-cutoff` is NOT `MAX`, tag it with an error: `contig not MAX` and skip the sequence (do not try to find primers)
  4. If the reference of the sample is `HIV1-CON-XX-Consensus-seed` tag that sequence with an error: `is V3 sequence`, skip the sequence (do not try to find primers), and write it to the `*primer_analysis.csv` file
  5. If there is an `X` in the middle of the sequence, tag that sequence with an error: `low internal read coverage`, skip the sequence (do not try to find primers), and write it to the `*primer_analysis.csv` file
  6. If there are ANY non-TCGA characters in the sequence, tag that sequence with an error: `contig sequence contained non-TCGA/gap`, skip the sequence (do not try to find primers), and write it to the `*primer_analysis.csv` file
  7. For each end (5' (fwd), 3' (rev)) of the sequence:
     1. If there are `X` characters found, try to remove them (if they are clustered) and if not possible to remove tag the fwd/rev end with a fwd/rev primer error: `low read coverage in primer region`, skip the fwd/rev end (do not try to find primers)
     2. If fwd/rev end has zero nucleotides found for primer, tag the fwd/rev end with a fwd/rev primer error: `primer was not found`, skip to the next end if any
     3. If the fwd/rev primer is deemed not valid, tag the fwd/rev end with a fwd/rev primer error: `primer failed secondary validation`, skip to the next end if any
  8. Write the sequence to the `*primer_analysis.csv` file regardless of tagged errors in any error column
  9. Load the `*primer_analysis.csv` files for both contigs and conseqs and for both of them apply the following filters in order:
     1. Remove all rows where either the `error`, `fwd_error`, or `rev_error` is tagged
     2. Remove the primers from the sequences (for hivseqinr)
     3. Remove rows where sample name appears twice (duplicates)
     4. Remove rows where the reference contains `unknown` or `reverse`
  10. Finally merge the filtered contigs and conseqs and write the final `*filtered.csv` file with conseqs taking precedence over contigs

### Creating `table_precursor.csv` from a fasta file
We got a request to create `table_precursor.csv`, containing only the gene regions trimmed out from the whole
sequence, based on a fasta file of sequences. I made a few modifications to allow for this. Specifically, follow 
these steps to do the same thing:
1. Run `gene_splicer.gene_splicer.py` on the input fasta. You may have to modify the sample name parsing [here][1] 
   to accommodate for the naming conventions in your input file. This should create folders for each sample within the
   specified results directory, each containing a file named `genes.fasta`, containing the sequences for all genes for 
   that sample.
2. To collate the information from the individual folders, run `generate_table_regions` from `gene_splicer.utils`,
   with your results folder as input argument. This is a modified version of `generate_table_precursor` and only 
   collates the gene region sequences. It should create a file named `table_precursor.csv` in your results folder!

[1]: https://github.com/cfe-lab/proviral/blob/46116021b4e70411cf6f0f0492e30292af83bd6e/gene_splicer/gene_splicer.py#L32