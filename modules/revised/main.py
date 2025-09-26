"""Healdette main module for antibody generation pipeline."""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List
from modules.revised.antibody_generator import AntibodyGenerator
from modules.revised.sequence_validator import SequenceValidator

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

def load_fusion_context(input_file: Path) -> Dict:
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

def save_results(results: Dict, output_file: Path):
    """Save generation results with formatting."""
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    logging.info(f"Results saved to {output_file}")

def generate_binders(fusion_context: Dict, num_candidates: int = 10) -> Dict:
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