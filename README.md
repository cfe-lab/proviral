
# CFE Proviral Pipeline

A comprehensive, Docker-enabled pipeline for proviral sequence
analysis. The CFE Proviral Pipeline automates the process of quality
filtering, primer detection and removal, sequence alignment, gene
splicing, and defect classification (using HIVSeqinR or CFEIntact) –
all built for ease of use by researchers with minimal programming
expertise.

Studying proviral genomes is critical for understanding viral
persistence, reactivation potential, and the barriers to curing
infections. The pipeline was developed to simplify this complex
analysis by:

- Ensuring high-quality sequence selection through automated filtering
  and control checks.
- Detecting laboratory-introduced primer sequences and removing them
  so that only genuine viral genomic data are analyzed.
- Extracting and aligning gene segments from the viral genome (based
  on the standard HIV HXB2 reference) and identifying key genetic
  defects.

In short, the pipeline supports researchers in gaining better insight
into the biology of viral reservoirs and the effects of antiviral
therapies.

---

## Features

- Primer Detection and Removal

- Quality Filtering
  - Applies multiple criteria to decide which sequences advance for analysis
  - Flags non-HIV, low read coverage, and sequences with internal ambiguities
  - Evaluates sequencing depth and coverage from MiCall's cascade output

- Sequence Alignment and Gene Splicing
  - Aligns sequences against HXB2 and extracts gene-level segments for further study

- Defect Detection and Prioritization
  - Integrates secondary modules (HIVSeqinR or CFEIntact) to assess genetic defects
  - Provides a "most severe defect" verdict for each sample in the final summary

- Comprehensive Output Summaries
  - Supports downstream visualization and statistical aggregation

- Docker-Based Deployment
  - Prepackaged Docker image reduces setup time and dependency conflicts

----

## Installation and Usage

Refer to [the documentation page](https://cfe-lab.github.io/proviral/introduction.html).

---

## Project Structure

The repository is organized as follows:

- `Dockerfile` – Defines build instructions for the main image (Ubuntu-based).
- `docs/` – Contains all documentation in markdown format (including installation, usage, workflow details, error codes, and troubleshooting).
- `cfeproviral/` – Main source code directory containing:
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
