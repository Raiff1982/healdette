import sys
from pathlib import Path
import os

# Add the test_output directory to Python path
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

from generate_binders_fixed import AntibodyGenerator

def main():
    # Initialize generator
    generator = AntibodyGenerator()
    
    # Test context with a real protein motif
    test_context = {
        "cleaned_sequence": "QVQLVQSGAEVKKPGASVKVSCKASGYTFT",  # Common VH framework motif
    }
    
    print("Initializing sequence generation test...")
    print(f"Target sequence: {test_context['cleaned_sequence']}")
    print("\nGenerating sequences...")
    
    # Generate 2 candidates to test
    results = generator.generate_binders(test_context, num_candidates=2)
    
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
        aa_counts = {aa: sequence.count(aa) for aa in set(sequence)}
        most_common = max(aa_counts, key=aa_counts.get)
        most_common_freq = aa_counts[most_common] / len(sequence)
        print(f"Most common AA: {most_common} ({most_common_freq:.1%})")
        
        # Check for homopolymer runs
        max_run = 1
        current_run = 1
        current_aa = sequence[0]
        for aa in sequence[1:]:
            if aa == current_aa:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
                current_aa = aa
        print(f"Longest homopolymer run: {max_run}")
        
        # Calculate hydrophobic content
        hydrophobic = "AILMFWYV"
        hydrophobic_count = sum(sequence.count(aa) for aa in hydrophobic)
        hydrophobic_fraction = hydrophobic_count / len(sequence)
        print(f"Hydrophobic content: {hydrophobic_fraction:.1%}")
        
        print("\nFull Sequence:")
        print(sequence)
        print("-" * 50)
    
    print(f"\nGeneration Statistics:")
    print(f"Success Rate: {results['stats']['success_rate']:.2%}")
    print(f"Total Attempts: {results['stats']['attempts']}")

if __name__ == "__main__":
    main()