from pathlib import Path
import sys
import os

# Add the project root directory to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from modules.generate_binders import AntibodyGenerator

def main():
    # Initialize generator with test data
    generator = AntibodyGenerator()
    
    # Override validation file path for testing
    generator.validation_set = []  # Start with empty validation set
    
    # Test context with a real protein motif
    test_context = {
        "cleaned_sequence": "QVQLVQSGAEVKKPGASVKVSCKASGYTFT",  # Common VH framework motif
    }
    
    print("Initializing sequence generation test...")
    print(f"Target sequence: {test_context['cleaned_sequence']}")
    print("\nGenerating sequences...")
    
    # Generate 3 candidates to test
    results = generator.generate_binders(test_context, num_candidates=3)
    
    # Print results
    print("\nGenerated Antibody Sequences:")
    print("=" * 50)
    
    for idx, binder in enumerate(results["generated_binders"], 1):
        print(f"\nCandidate {idx}:")
        print(f"Template VH: {binder['template_vh']}")
        print(f"Template VL: {binder['template_vl']}")
        print(f"Validation Score: {binder['validation_score']:.2f}")
        print("\nHeavy Chain CDRs:")
        for i, cdr in enumerate(binder['heavy_cdrs'], 1):
            print(f"CDR-H{i}: {cdr}")
        print("\nLight Chain CDRs:")
        for i, cdr in enumerate(binder['light_cdrs'], 1):
            print(f"CDR-L{i}: {cdr}")
            
        # Print additional analysis
        sequence = binder['sequence']
        print("\nSequence Analysis:")
        print(f"Length: {len(sequence)} AA")
        aa_counts = {aa: sequence.count(aa) / len(sequence) for aa in set(sequence)}
        print(f"Max AA frequency: {max(aa_counts.values()):.2%}")
        print(f"Most common AA: {max(aa_counts, key=aa_counts.get)}")
        
        # Check for homopolymer runs
        max_run = 1
        current_run = 1
        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
        print(f"Longest homopolymer run: {max_run}")
        
        print("\nFull Sequence:")
        print(binder['sequence'])
        print("-" * 50)
    
    print(f"\nGeneration Statistics:")
    print(f"Success Rate: {results['stats']['success_rate']:.2%}")
    print(f"Total Attempts: {results['stats']['attempts']}")

if __name__ == "__main__":
    main()