---
title: Sample Analysis Process
---

This document explains how the "sample" entrypoint (implemented in `sample.py`) orchestrates the analysis --- from preprocessing the raw sequences to producing the final output files. It describes how the pipeline reads input data, applies primer detection and gene splicing, and generates summary files ready for downstream use.

---

## Overview

When you run the command--line entrypoint (for example, using the Docker image with the command `cfeproviral sample ...`), the pipeline follows these major steps:

1. **Command-Line & Environment Setup**
2. **Primer Finding and Filtering**
3. **Gene Splicing (Extracting Viral Genes)**
4. **Generating an Upload-Ready Table**
5. **Building the Proviral Landscape File**
6. **Final Output File Consolidation**

Below, each step is described in detail.

---

## 1. Command-Line Parsing & Setup

- **Input Files:**
  The script accepts several positional parameters:
  - `sample_info.csv` --- Contains metadata (run name, sample name, etc.).
  - `contigs.csv` and `conseqs.csv` --- Hold raw sequences from assembled contigs and consensus sequences.
  - `cascade.csv` --- Produced by MiCall, this file provides a read count (remap) used to decide if a sequence should be analyzed.

- **Output Files:**
  You must supply paths for several outputs:
  - An overall **Outcome Summary** CSV.
  - Two primer analysis CSVs (one for consensus, another for contigs).
  - A **Table Precursor** CSV (ready for downstream analysis).
  - A **Proviral Landscape** CSV (for plotting gene-level defects).
  - An archive (tar) of **Detailed Results** produced by the secondary analysis (HIVSeqinR or CFEIntact).

- **Additional Options:**
  Options such as the probe size (i.e. the number of nucleotides examined at each sequence end) and flags to choose the secondary analysis backend are also processed.

- **Metadata Parsing:**
  The pipeline reads `sample_info.csv` to extract the run and sample names. These values are used throughout the process and to label output files.

---

## 2. Primer Finding and Filtering

- **Primary Function:**
  The primer finder module (invoked via `primer_finder.run()`) handles the following:

  - **Grouping Sequences:**
    Using the `cascade.csv`, the pipeline groups data by sample and checks that each sample has enough reads to be analyzed.

  - **Searching for Primers:**
    It examines a fixed--length "probe" region from the ends of each sequence (from the consensus or contig CSVs) to look for the correct primer sequences (forward and reverse).
    - If the expected primer isn't found, the pipeline may try the reverse complement of the sequence.
    - Error conditions (e.g., low coverage, non--TCGA characters, or missing primers) are recorded with specific error codes.

  - **Output:**
    The results are written to a primer analysis CSV file. A merged "joined" CSV file is also produced that combines information from the contig and consensus analyses.
    Finally, the primers are "stripped" from the sequences and one or more FASTA files are generated. These FASTA files contain the final, cleaned sequence (often with synthetic primer sequences added back if required by the secondary analysis).

---

## 3. Gene Splicing (Extracting Viral Genes)

- **Alignment:**
  For each FASTA file produced by the primer finder, the gene splicer module is invoked via `gene_splicer.run()`. This module:

  - Uses an external aligner (called `minimap2`) to align the processed sample sequence against a modified version of the HXB2 reference genome.

- **Extracting Genes:**
  Based on the alignment and using predefined gene annotations:
  - It "splices" the viral genes by mapping the alignment coordinates to annotated gene regions.
  - The extracted gene sequences are then written to a file (usually named `genes.fasta` in the sample's folder).

This step makes it possible to later perform gene--level analysis (for example, checking genomic defects in gag, pol, env, etc.).

---

## 4. Generating an Upload-Ready Table

- **Merging Results:**
  Once primer filtering and gene splicing are complete, the pipeline combines the final processed sequences with the defect classification (obtained from downstream analysis tools such as HIVSeqinR or CFEIntact).

- **Table Precursor:**
  The utility function `generate_table_precursor()` produces a file named `table_precursor.csv`. This file aggregates:

  - The cleaned, primer--stripped sequence.
  - Meta-information about the sequence (such as its type and length).
  - Defect calls and gene--level regions (as extracted during gene splicing).

This table is designed to be used for subsequent analysis steps or uploaded directly into visualization tools.

---

## 5. Building the Proviral Landscape File

- **Landscape Generation:**
  The landscapes module takes results from the secondary analysis backend and, by reading BLAST outputs or other alignment reports, creates a CSV file (`proviral_landscape.csv`). This file includes:
  - Genomic coordinates (relative to the reference genome) for the different fragments of the proviral sequence.
  - Flags indicating the presence of defects (or inversions) in those fragments.

- **Purpose:**
  The landscape file is intended for visualizing the structure of the provirus. Researchers can use it to understand the distribution of defects across the viral genome.

---

## 6. Final Output File Consolidation

- **File Copying:**
  After processing is complete and all intermediate files exist in a temporary "scratch" directory, the script moves or renames the outputs into their final user-designated locations. These include:

  - Outcome summary (overall sample report).
  - Consensus and contig primer analysis files.
  - The table precursor for downstream upload.
  - The proviral landscape CSV.
  - A tar archive with detailed results from the secondary analysis.

- **Completion:**
  Once these files are in place, the pipeline completes its run and exits. Users can then inspect the outputs for downstream analysis or troubleshooting.

---

## Summary

In short, the "sample" entrypoint provided by the pipeline performs the following sequence of actions:

- **Parse and validate inputs** (raw sequences, metadata, MiCall cascade data).
- **Find, validate, and strip primers** from each sample sequence while logging any error conditions.
- **Align sequences to the HXB2 reference** and extract the viral gene segments.
- **Merge and summarize the results** into an overall outcome summary and a table precursor file.
- **Generate a proviral landscape file** that visualizes genomic defects.
- **Consolidate all outputs** into user--specified final files for use in further analysis.

By encapsulating this entire process in a Docker container, the pipeline minimizies software--dependency issues and makes the sample analysis accessible to users with minimal programming expertise.
