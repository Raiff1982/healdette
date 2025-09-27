"""
Test script for generating and validating antibody sequences.
"""
from modules.sequence_generator import SequenceGenerator
from modules.simple_validator import SimpleValidator
import json

def format_validation_results(sequence, results):
    """Format validation results for display."""
    output = [
        f"Sequence: {sequence}",
        f"Valid: {results['valid']}",
        f"\nMetrics:",
        f"  Aromatic content: {results['metrics']['aromatic_content']:.1f}%",
        f"  Hydrophobic content: {results['metrics']['hydrophobic_content']:.1f}%",
        f"  Net charge: {results['metrics']['net_charge']}",
        f"\nPopulation Scores:"
    ]
    
    for pop, score in results['population_scores'].items():
        output.append(f"  {pop}: {score['score']:.2f} (weight: {score['weight']:.3f})")
    
    if results['warnings']:
        output.extend(["\nWarnings:", *[f"  - {w}" for w in results['warnings']]])
    
    return "\n".join(output)

def main():
    config_file = "config_test.json"
    
    # Initialize generator and validator
    generator = SequenceGenerator(config_file)
    validator = SimpleValidator(config_file)
    
    # Generate sequences
    print("Generating antibody sequences...")
    sequences = generator.generate_sequences(num_sequences=5)
    
    # Validate each sequence
    print("\nValidating sequences against multi-ethnic parameters...")
    print("-" * 80)
    
    valid_sequences = []
    for i, sequence in enumerate(sequences, 1):
        print(f"\nSequence {i}:")
        results = validator.validate_sequence(sequence)
        print(format_validation_results(sequence, results))
        print("-" * 80)
        
        if results['valid']:
            valid_sequences.append({
                'sequence': sequence,
                'validation': results
            })
    
    print(f"\nSummary:")
    print(f"Total sequences generated: {len(sequences)}")
    print(f"Valid sequences: {len(valid_sequences)}")
    
    # Save valid sequences to file
    if valid_sequences:
        output_file = "valid_sequences.json"
        with open(output_file, 'w') as f:
            json.dump(valid_sequences, f, indent=2)
        print(f"\nValid sequences saved to {output_file}")

if __name__ == "__main__":
    main()