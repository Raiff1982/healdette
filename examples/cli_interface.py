"""
CLI example for Healdette antibody sequence generation pipeline.
"""
import argparse
import json
from pathlib import Path
from modules.weighted_validator import WeightedSequenceValidator
from modules.config_validator import ConfigValidator
from modules.sequence_generator import SequenceGenerator

def parse_args():
    parser = argparse.ArgumentParser(description='Healdette CLI Interface')
    parser.add_argument('config', type=str, help='Path to configuration JSON file')
    parser.add_argument('output', type=str, help='Path to output JSON file')
    parser.add_argument('--num-sequences', type=int, default=5,
                       help='Number of sequences to generate (default: 5)')
    parser.add_argument('--min-score', type=float, default=0.7,
                       help='Minimum validation score to include sequence (default: 0.7)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load and validate configuration
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: Configuration file {args.config} not found")
        return 1
    
    config_validator = ConfigValidator()
    validation_result = config_validator.validate_file(args.config)
    
    if not validation_result['valid']:
        print("Error: Invalid configuration")
        for error in validation_result['errors']:
            print(f"  - {error}")
        return 1
    
    # Generate sequences
    generator = SequenceGenerator(args.config)
    sequences = generator.generate_sequences(num_sequences=args.num_sequences)
    
    # Validate sequences
    validator = WeightedSequenceValidator(args.config)
    results = []
    
    for seq in sequences:
        validation = validator.validate_sequence(seq)
        if validation['valid']:
            # Calculate overall score
            pop_scores = validation['population_scores']
            weighted_score = sum(s['score'] * s['weight'] for s in pop_scores.values())
            
            if weighted_score >= args.min_score:
                results.append({
                    'sequence': seq,
                    'validation': validation,
                    'weighted_score': weighted_score
                })
    
    # Write results to output file
    with open(args.output, 'w') as f:
        json.dump({
            'sequences': results,
            'total_generated': args.num_sequences,
            'total_valid': len(results),
            'min_score': args.min_score
        }, f, indent=2)
    
    print(f"\nGenerated {len(results)} valid sequences")
    print(f"Results written to {args.output}")
    return 0

if __name__ == '__main__':
    exit(main())