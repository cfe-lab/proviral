#!/bin/bash

# I think R needs this directory?
mkdir -p /usr/share/doc/R-3.6.0/html/

local({r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})

message('======== INSTALLING BIOSTRINGS ========')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

message('======== INSTALLING MUSCLE ========')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("muscle")