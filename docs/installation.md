---
title: Installation
---

To streamline the installation process and avoid potential software conflicts, we utilize Docker. 
Docker allows us to package all necessary components inside a container, ensuring consistency across different environments.

## Install Docker

- **Windows / macOS**:
  1. Visit the [Docker Desktop download page](https://www.docker.com/products/docker-desktop/) and download the installer for your operating system.
  2. Follow the installation instructions provided on Docker's website.
  3. Once installed, launch Docker Desktop and ensure it is running.

- **Linux**:
  1. Follow the [Docker Linux installation guide](https://docs.docker.com/engine/install/) for your distribution (e.g., Ubuntu, CentOS).
  2. Start the Docker service using:
     ```bash
     sudo systemctl start docker
     ```
  3. Optionally, enable Docker to start at boot:
     ```bash
     sudo systemctl enable docker
     ```

## Pull the Docker Image for Proviral Analysis

With Docker installed and running, you can now download the pre-configured Docker image that contains the proviral sequence analysis pipeline:

- Open your terminal or command prompt and execute the following command:

  ```bash
  docker pull cfelab/proviral
  ```

  This command downloads the latest version of the Docker image from the repository.
  To test the image, execute:
  
  ```bash
  docker run --rm cfelab/proviral --version
  ```
  
  You should see a version number like this:
  
  ```
  Proviral Pipeline 2.4.0
  ```

---

Next: [data preparation](data_prep.html).
