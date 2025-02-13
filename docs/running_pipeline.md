---
title: Running the Pipeline
---

The proviral pipeline is designed to analyze proviral sequences through several stages --- from primer detection and alignment to detailed error reporting and result summarization. In this document we assume you have already installed Docker and have prepared your input data files. In this chapter you will learn how to run the pipeline using the provided Docker image as well as the sample entrypoint for a complete end--to--end example.

---

## Assumptions and Prerequisites

Before you start, please ensure that:

- Docker is installed and running on your system.
- Your current working directory contains the required input files:
  - **sample_info.csv** (a CSV containing sample metadata such as run_name and sample identifiers)
  - **contigs.csv** (assembled contig sequences)
  - **conseqs.csv** (consensus sequences)
  - **cascade.csv** (MiCall output indicating mapped read counts)

For more details on data preparation and installation, see the [Installation](installation.html) and [Data Preparation](data_prep.html) sections.

---

## Launching the Pipeline Using Docker

The pipeline image available on Docker Hub is preconfigured to run all necessary steps. In this example, we demonstrate how to run the pipeline through the `sample` entrypoint. This entrypoint coordinates the analysis by calling several modules (primer detection, gene splicing, landscape generation, etc.) and ultimately writes a set of output files for review.

To run the pipeline using the `sample` entrypoint, open your terminal, navigate to your working directory (which contains your input files), and use the following command:

```bash
docker run --rm -v .:/w cfelab/proviral sample sample_info.csv contigs.csv conseqs.csv cascade.csv outcome_summary.csv conseqs_primers.csv contigs_primers.csv table_precursor.csv proviral_landscape.csv detailed_results.tar --cfeintact
```

### Command Breakdown

- **`docker run --rm`**
  Runs the container and automatically removes it after execution.

- **`-v .:/w`**
  Mounts your current directory (which holds both input and eventual output files) into the container’s `/w` folder.

- **`cfelab/proviral`**
  Specifies the Docker image containing the proviral pipeline.

- **`sample`**
  Instructs the container to use the `sample` entrypoint (the *sample.py* module). This entrypoint is ideal for processing a single-sample run.

- **Positional Parameters:**
  `sample_info.csv contigs.csv conseqs.csv cascade.csv outcome_summary.csv conseqs_primers.csv contigs_primers.csv table_precursor.csv proviral_landscape.csv detailed_results.tar`
  These parameters denote, in order:
  - **sample_info.csv** – An input metadata file
  - **contigs.csv** and **conseqs.csv** – Input sequence files
  - **cascade.csv** – The cascade CSV from MiCall
  - **Output Files:**
    - **outcome_summary.csv** – Overall summary of the analysis
    - **conseqs_primers.csv** – Primer analysis for consensus sequences
    - **contigs_primers.csv** – Primer analysis for contigs
    - **table_precursor.csv** – Data ready for downstream upload
    - **proviral_landscape.csv** – Data for generating proviral landscape plots
    - **detailed_results.tar** – Archive of detailed results produced by HIVSeqinR (or CFEIntact)

- **`--hivseqinr`**
  Specifies that the HIVSeqinR backend should be used for downstream analysis. (You may alternatively specify `--cfeintact` if that is preferred.)

---

## Troubleshooting

In case you encounter errors:

- Review the terminal output and any generated log files to identify issues (common errors include missing or misformatted input files).
- Verify that Docker is properly mounting your current directory (using the `-v` flag).
- Revisit [Installation](installation.html) and [Data Preparation](data_prep.html) if input file formats or configurations are uncertain.
- Run the pipeline with the `--help` option inside Docker (e.g., `docker run --rm cfelab/proviral --help`) to review available options and usage information.

---

Next: [Viewing and Interpreting Outputs](interpretation.html).
