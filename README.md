
# Healdette: Antibody Sequence Generation Pipeline

A computational pipeline for generating and validating antibody sequences using machine learning and biophysical analysis. The pipeline integrates ProtGPT2 for sequence generation with BioPython for structural analysis and includes population-specific HLA frequency data for immunogenicity assessment.

## Features

- Antibody sequence generation using ProtGPT2 with template-based constraints
- Sequence validation against known therapeutic antibodies
- Population-specific immunogenicity assessment using HLA frequency data
- Biophysical property analysis using BioPython
- Structured output in JSON format with detailed analysis results

## Requirements

- Python 3.8 or higher
- CUDA-capable GPU (recommended for ProtGPT2)
- Required Python packages listed in `requirements.txt`

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Raiff1982/healdette.git
cd healdette
```

2. Create and activate a virtual environment:
```bash
python -m venv .venv
# On Windows:
.venv\Scripts\activate
# On Unix/MacOS:
source .venv/bin/activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

1. Configure input parameters in a JSON file:
```json
{
    "target_sequence": "EVQLVESGGGLVQPGGSLRLSCAASGFTFS",
    "population": "EUR",
    "num_sequences": 10
}
```

2. Run the pipeline:
```bash
python main.py --config input_config.json
```

## Output Files

The pipeline generates two types of output files in the `output` directory:

1. Detailed JSON output (`antibody_designs_{timestamp}.json`):
   - Generated antibody sequences with framework and CDR regions
   - Biophysical properties (hydrophobicity, charge, stability)
   - Population-specific immunogenicity scores
   - Validation results against therapeutic antibodies

2. Summary report (`antibody_summary_{timestamp}.txt`):
   - Key metrics for each generated sequence
   - Population coverage statistics
   - Validation summary

## Reproducibility

To reproduce the results:

1. Use the same random seed for ProtGPT2:
```python
import torch
torch.manual_seed(42)
```

2. Ensure consistent data sources:
   - HLA frequency data: NetMHCpan 4.1 database
   - Therapeutic antibody dataset: THERAb database v2.0
   - Framework templates: IMGT database

3. Run validation tests:
```bash
python -m unittest discover tests
```

## License

MIT License. See LICENSE file for details.

## Citation

If you use this software in your research, please cite:
```
Harrison, J. (2025). Healdette: A Population-Aware Antibody Design Pipeline.
GitHub repository: https://github.com/Raiff1982/healdette
```

## Author

Jonathan Harrison (Raiff1982)
