"""
Basic usage example for Healdette antibody sequence generation pipeline.
"""
from modules.weighted_validator import WeightedSequenceValidator
from modules.config_validator import ConfigValidator
from modules.sequence_generator import SequenceGenerator

def main():
    # Load and validate configuration
    config_path = "european_populations_config.json"
    config_validator = ConfigValidator()
    
    if config_validator.validate_file(config_path)['valid']:
        # Initialize sequence generator
        generator = SequenceGenerator(config_path)
        
        # Generate candidate sequences
        sequences = generator.generate_sequences(num_sequences=5)
        
        # Create validator
        validator = WeightedSequenceValidator(config_path)
        
        # Validate and score each sequence
        for i, seq in enumerate(sequences, 1):
            results = validator.validate_sequence(seq)
            
            print(f"\nSequence {i}:")
            print(f"Sequence: {seq}")
            print(f"Valid: {results['valid']}")
            
            if results['valid']:
                print("\nPopulation Scores:")
                for pop, score in results['population_scores'].items():
                    print(f"  {pop}: {score['score']:.2f} (weight: {score['weight']:.3f})")
                
                print("\nMetrics:")
                metrics = results['metrics']
                print(f"  Aromatic content: {metrics['aromatic_content']:.1f}%")
                print(f"  Hydrophobic content: {metrics['hydrophobic_content']:.1f}%")
                
                binding_scores = metrics['binding_motifs']['scores']
                print("\nBinding Motif Scores:")
                for pop, score in binding_scores.items():
                    print(f"  {pop}: {score['score']:.2f} (weighted: {score['weighted_score']:.3f})")
                print(f"  Total binding score: {metrics['binding_motifs']['total_score']:.3f}")
            
            if results['warnings']:
                print("\nWarnings:")
                for warning in results['warnings']:
                    print(f"  - {warning}")

if __name__ == "__main__":
    main()