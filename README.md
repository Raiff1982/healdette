
# Healdette: Antibody Sequence Generation Pipeline

A computational pipeline for generating and validating antibody sequences using machine learning and biophysical analysis. The pipeline integrates ProtGPT2 for sequence generation with BioPython for structural analysis and includes multi-ethnic HLA frequency data for immunogenicity assessment, with optimizations for various population-specific binding motifs.

## Features

- Antibody sequence generation using ProtGPT2 with template-based constraints
- Multi-ethnic binding motif optimization:
  - Celtic: WY, RF, KW, YF pairs
  - Asian: WH, RY, KF, HF pairs
  - Mediterranean: FY, RW, KY, WF pairs
- Population-specific sequence validation parameters:
  Celtic:
    - Aromatic content: 15-27%
    - Hydrophobic content: 35-45%
    - Net charge: +5 to +15
  Asian:
    - Aromatic content: 12-25%
    - Hydrophobic content: 30-40%
    - Net charge: +3 to +12
  Mediterranean:
    - Aromatic content: 18-30%
    - Hydrophobic content: 32-42%
    - Net charge: +4 to +14
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

## Multi-Ethnic Configuration

Healdette now supports ancestry-weighted validation for multiple ethnic populations. The system uses:
1. Population-specific binding motifs and parameters
2. Ancestry weights from genetic analysis
3. HLA frequency data for immunogenicity assessment

### Configuration Structure

Configuration files follow this structure:
```json
{
    "global_params": {
        "sequence_length": {
            "min": 40,
            "max": 70
        },
        "structural_params": {
            "helix_propensity": {
                "min": 20,
                "max": 50
            },
            "sheet_propensity": {
                "min": 10,
                "max": 40
            }
        },
        "homopolymer_threshold": 4
    },
    "populations": {
        "french_german": {
            "ancestry_weight": 0.298,
            "binding_motifs": ["WY", "RF", "KH", "YF"],
            "biophysical_params": {
                "aromatic_content": {
                    "min": 16,
                    "max": 28
                },
                "hydrophobic_content": {
                    "min": 33,
                    "max": 43
                },
                "net_charge": {
                    "min": 4,
                    "max": 13
                }
            },
            "hla_frequencies": {
                "hla_a": {},
                "hla_b": {},
                "hla_c": {}
            }
        }
    }
}
```

### Ancestry-Weighted Validation

The validation system considers:
1. **Ancestry Weights**: Each population's contribution is weighted by ancestry percentage
2. **Blended Parameters**: Biophysical parameters are blended based on ancestry weights
3. **Multiple Binding Motifs**: Scores binding motifs from all relevant populations
4. **HLA Compatibility**: Considers population-specific HLA frequencies

### Population-Specific Parameters

Each population can define:
- **Binding Motifs**: Amino acid pairs crucial for binding
- **Biophysical Parameters**:
  - Aromatic content ranges
  - Hydrophobic content ranges
  - Net charge requirements
- **HLA Frequencies**: Population-specific HLA allele distributions

## Usage

1. Create a configuration file following the schema (see `examples/` directory):
```json
{
    "global_params": {
        "sequence_length": {
            "min": 40,
            "max": 70
        }
    },
    "populations": {
        "french_german": {
            "ancestry_weight": 0.298,
            "binding_motifs": ["WY", "RF", "KH", "YF"],
            "biophysical_params": {
                "aromatic_content": {
                    "min": 16,
                    "max": 28
                }
            }
        },
        "finnish": {
            "ancestry_weight": 0.057,
            "binding_motifs": ["WH", "RF", "KY", "FF"],
            "biophysical_params": {
                "aromatic_content": {
                    "min": 14,
                    "max": 26
                }
            }
        }
    }
}
```

2. Validate sequences using the weighted validator:
```python
from modules.weighted_validator import WeightedSequenceValidator
from modules.config_validator import ConfigValidator

# Load and validate configuration
config_validator = ConfigValidator()
config = "path/to/config.json"
if config_validator.validate_file(config)['valid']:
    # Create validator with ancestry-weighted parameters
    validator = WeightedSequenceValidator(sequence, config)
    
    # Get detailed validation results
    results = validator.validate_sequence()
    
    # Check population-specific scores
    pop_scores = results['population_scores']
    for pop, score in pop_scores.items():
        print(f"{pop}: {score['score']} (weight: {score['weight']})")
```

3. Run the pipeline with multi-ethnic configuration:
```bash
python main.py config.json output.json --num-candidates 15
```

### Example Configurations

Complete example configurations are available in the `examples/` directory:
- `european_populations_config.json`: Configuration for European population clusters
- `multi_ethnic_config.json`: General multi-ethnic configuration template
- `celtic_test_input.json`: Celtic-specific test configuration

### Understanding Validation Results

The weighted validator provides detailed results:
```json
{
    "valid": true,
    "warnings": [],
    "metrics": {
        "aromatic_content": 22.5,
        "hydrophobic_content": 38.2,
        "binding_motifs": {
            "scores": {
                "french_german": {"score": 0.75, "weighted_score": 0.223},
                "finnish": {"score": 0.5, "weighted_score": 0.029}
            },
            "total_score": 0.252
        }
    },
    "population_scores": {
        "french_german": {
            "score": 0.8,
            "weight": 0.298
        },
        "finnish": {
            "score": 0.6,
            "weight": 0.057
        }
    }
}
    ],
    "num_sequences": 10,
    "global_validation_params": {
        "min_sequence_length": 40,
        "max_sequence_length": 70,
        "allow_homopolymers": false,
        "structure_requirements": {
            "helix_propensity": {
                "min": 0.2,
                "max": 0.5
            },
            "sheet_propensity": {
                "min": 0.1,
                "max": 0.4
            }
        }
    }
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
   - Celtic binding motif analysis
   - Biophysical properties (hydrophobicity, charge, stability)
   - Aromatic content and distribution
   - Population-specific immunogenicity scores
   - Validation results against therapeutic antibodies

2. Summary report (`antibody_summary_{timestamp}.txt`):
   - Key metrics for each generated sequence
   - Celtic motif occurrence statistics
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
   - Celtic binding motif templates: Custom database

3. Run validation tests:
```bash
python -m unittest discover tests
```

## License

MIT License. See LICENSE file for details.

## Citation

If you use this software in your research, please cite:
```bibtex
@software{healdette2025,
  title = {Healdette: Celtic-Optimized Antibody Generation Pipeline},
  author = {Raiff, et al.},
  year = {2025},
  version = {1.0.0},
  url = {https://github.com/Raiff1982/healdette}
}
```
Harrison, J. (2025). Healdette: A Population-Aware Antibody Design Pipeline.
GitHub repository: https://github.com/Raiff1982/healdette
```

## Author

Jonathan Harrison (Raiff1982)
