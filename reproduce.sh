#!/bin/bash
set -e

# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate

# Install minimal dependencies
pip install -r requirements.txt

# Run validation with deterministic mode
python run_pipeline.py --deterministic

# Generate visualization
python visualize_results.py

# Verify checksums
sha256sum -c checksums.sha256

echo "Reproduction completed successfully!"