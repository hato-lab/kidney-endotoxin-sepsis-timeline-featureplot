#!/usr/bin/env bash
set -euo pipefail

# This script downloads the necessary R data object file from GEO.

# Download and gunzip the R data object from GEO.
mkdir -p data
wget -O /dev/stdout ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151658/suppl/GSE151658%5Fintegrated%2E0h%2E1h%5F4h%5F16h%5F27h%5F36h%5F48h%2Erds%2Egz | gunzip -c > data/integrated.0h.1h_4h_16h_27h_36h_48h.rds
