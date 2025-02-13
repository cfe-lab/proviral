
# CFE Proviral Pipeline

A comprehensive, Docker-enabled pipeline for proviral sequence analysis. The CFE Proviral Pipeline automates the process of quality filtering, primer detection and removal, sequence alignment, gene splicing, and defect classification (using HIVSeqinR or CFEIntact) – all built for ease of use by researchers with minimal programming expertise.

Proviral genomes — in particular, those integrated into host DNA during HIV infections — are a key subject of study in viral persistence research. Analyzing these sequences can reveal both intact and defective viral genomes, with implications for HIV reactivation and cure strategies.

---

## Features

- Primer Detection and Removal
- Quality Filtering
  - Evaluates sequencing depth and coverage from MiCall's cascade output
  - Flags non-HIV, low read coverage, and sequences with internal ambiguities
  - Applies multiple criteria to decide which sequences advance for analysis

- Sequence Alignment and Gene Splicing
  - Extracts gene-level segments for further analysis

- Defect Detection and Prioritization
  - Integrates secondary modules (HIVSeqinR or CFEIntact) to assess genetic defects
  - Provides a "most severe defect" verdict for each sample in the final summary

- Comprehensive Output Summaries
  - Generates files such as `outcome_summary.csv`, `table_precursor.csv`, and `proviral_landscape.csv`
  - Supports downstream visualization and statistical aggregation

- Docker-Based Deployment
  - Prepackaged Docker image reduces setup time and dependency conflicts
  - Consistent runtime across diverse computing environments

----

## Installation and Usage

Refer to [the documentation page](https://cfe-lab.github.io/proviral/introduction.html).

---

## Project Structure

The repository is organized as follows:

- `Dockerfile` – Defines build instructions for the main image (Ubuntu-based).
- `docs/` – Contains all documentation in markdown format (including installation, usage, workflow details, error codes, and troubleshooting).
- `gene_splicer/` – Main source code directory containing:
  - `primer_finder.py`: Primer detection, filtering, and error logging.
  - `gene_splicer.py`: Gene splicing and alignment with external tools.
  - `main.py`: Entry point (cfeproviral) and command-line interface orchestration.
  - Additional modules, helpers, and utilities to support functions such as statistics, failure summary, and landscape generation.
- `setup.py` – Package installation script with dependency management and console-script entry point.

---

## Contributing

Contributions are welcome! Please review the following guidelines before submitting pull requests:

- Report bugs and feature requests via the [GitHub Issues tracker](https://github.com/cfe-lab/proviral/issues/new).
- Follow the established coding style and add tests where appropriate.
- Update documentation (both in docs/ and the README) for any changes that affect usage or behavior.
- See [the contributing file](docs/contributing.md) for a detailed guide on how to contribute.

---

## References and Acknowledgments

- Thanks to the developers behind HIVSeqinR and CFEIntact – their tools are integrated as backends for defect analysis.
- Acknowledgment to the MiCall project for alignment modules and related utilities.
- For further reading on proviral sequencing and HIV research methodologies, refer to the [pipeline's documentation](https://cfe-lab.github.io/proviral).
