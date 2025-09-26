"""Main module for antibody generation pipeline."""

import argparse
import logging
from pathlib import Path
import json
from antibody_generator import AntibodyGenerator
from sequence_validator import SequenceValidator

def setup_logging():
    """Configure logging for the application."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('healdette.log'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def load_fusion_context(input_file: Path) -> dict:
    """Load and validate fusion protein context."""
    try:
        with open(input_file) as f:
            context = json.load(f)
        
        required_fields = ['cleaned_sequence', 'fusion_points']
        missing_fields = [field for field in required_fields if field not in context]
        if missing_fields:
            raise ValueError(f"Missing required fields: {', '.join(missing_fields)}")
        
        return context
    except (json.JSONDecodeError, FileNotFoundError) as e:
        raise ValueError(f"Error loading fusion context from {input_file}: {str(e)}")

def save_results(results: dict, output_file: Path):
    """Save generation results with formatting."""
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    logging.info(f"Results saved to {output_file}")

def generate_binders(fusion_context: dict, num_candidates: int = 10) -> dict:
    """Generate and validate antibody binders."""
    generator = AntibodyGenerator()
    validator = SequenceValidator()
    
    # Generate initial candidates
    results = generator.generate_binders(fusion_context, num_candidates)
    
    # Additional validation of complete sequences
    for binder in results['generated_binders']:
        sequence = binder['sequence']
        validation_result = validator.analyze_sequence(sequence)
        binder['detailed_validation'] = validation_result
    
    # Personalize binders with HLA compatibility
    patient_data = {
        "immune_profile": [
            "A*01:01", "A*03:01",  # Celtic markers
            "B*07:02", "B*08:01",  # British Isles common
            "C*07:01", "C*04:01"   # Celtic-specific
        ],
        "metabolic_rate": 1.0,  # Standard for Northern European populations
        "prior_exposure": ["SARS-CoV-2", "Influenza-A", "EBV", "CMV"],
        "ancestry_profile": ["Celtic_British"],
        "population_weights": {
            "irish": 0.42,
            "scottish": 0.28,
            "english": 0.18,
            "germanic": 0.08,
            "welsh": 0.04
        }
    }
    
    if results['generated_binders']:
        from .personalize_binders import personalize_binders
        results['generated_binders'] = personalize_binders(
            results['generated_binders'],
            patient_data
        )
    
    return results

def main():
    """Main execution flow."""
    parser = argparse.ArgumentParser(description='Generate antibody binders for fusion proteins')
    parser.add_argument('input_file', type=Path, help='Input JSON file with fusion context')
    parser.add_argument('output_file', type=Path, help='Output JSON file for results')
    parser.add_argument('--num-candidates', type=int, default=10, help='Number of candidates to generate')
    args = parser.parse_args()
    
    logger = setup_logging()
    logger.info("Starting Healdette antibody generation pipeline")
    
    try:
        # Load and validate input
        fusion_context = load_fusion_context(args.input_file)
        logger.info(f"Loaded fusion context from {args.input_file}")
        
        # Generate binders
        results = generate_binders(fusion_context, args.num_candidates)
        logger.info(f"Generated {len(results['generated_binders'])} binders")
        
        # Save results
        save_results(results, args.output_file)
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        raise
    
    logger.info("Pipeline completed successfully")

if __name__ == "__main__":
    main()