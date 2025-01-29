
# Introduction

In the study of viral sequences, understanding the proviral genome -- the viral genome integrated into the host DNA -- is crucial. This guide is designed to help biology researchers without programming expertise analyze proviral sequences efficiently using a dedicated bioinformatics pipeline. By leveraging a Docker container, we simplify the technical setup, allowing you to focus on the biological insights from your data.

---

### Background

Proviral sequences are fragments of viral DNA integrated into the host's genome. Analyzing these sequences helps researchers understand viral integration, replication fidelity, and structure of proviral reservoir. Our pipeline facilitates this analysis by automating processes such as primer identification, sequence alignment, and error detection.

---

### Pipeline Features

The pipeline offers the following core features:

- **Primer Analysis**: Identifies primer sequences crucial for polymerase chain reactions.
- **Sequence Alignment and Filtering**: Uses tools like `minimap2` to align sequences and applies filters to ensure data quality.
- **Error Detection**: Categorizes sequence errors, such as non-proviral sequences or low coverage.
- **Result Summarization**: Generates comprehensive reports to aid in data interpretation and decision-making.

---

## Setup and Environment

To make our pipeline accessible and easy to use, we have encapsulated it within a Docker container. This ensures all dependencies and configurations are pre-installed, minimizing compatibility issues across different systems.

- **Docker**: A platform to develop, ship, and run applications in isolated environments. Ensure that Docker is installed and running on your computer.

By following the rest of this guide, you will be up and running with proviral sequence analysis quickly, even with minimal technical background.

---

Next: [installation](/installation.html).
