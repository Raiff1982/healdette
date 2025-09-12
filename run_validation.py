"""
Run sequence validation on the generated antibody designs.
"""

from modules.validate_sequences import validate_binder_set
import json
from datetime import datetime

def analyze_validation_results(results):
    """Analyze and print validation results summary."""
    binders = results["validated_binders"]
    
    print("\nValidation Summary")
    print("-" * 50)
    
    # Disorder analysis
    disorder_scores = [b["validation"]["disorder"] for b in binders]
    print(f"\nDisorder Analysis:")
    print(f"Average disorder: {sum(disorder_scores)/len(disorder_scores):.3f}")
    print(f"Range: {min(disorder_scores):.3f} - {max(disorder_scores):.3f}")
    
    # Signal peptides
    signal_count = sum(1 for b in binders if b["validation"]["signal_peptide"]["has_signal"])
    print(f"\nSignal Peptide Detection:")
    print(f"Sequences with potential signal peptides: {signal_count}/{len(binders)}")
    
    # Cysteine analysis
    cys_counts = [b["validation"]["cysteines"]["count"] for b in binders]
    paired_cys = sum(1 for b in binders if b["validation"]["cysteines"]["paired"])
    print(f"\nCysteine Analysis:")
    print(f"Sequences with paired cysteines: {paired_cys}/{len(binders)}")
    print(f"Cysteine count range: {min(cys_counts)} - {max(cys_counts)}")
    
    # Glycosylation sites
    glyco_sites = [len(b["validation"]["glycosylation"]) for b in binders]
    print(f"\nGlycosylation Sites:")
    print(f"Total potential sites: {sum(glyco_sites)}")
    print(f"Average sites per sequence: {sum(glyco_sites)/len(glyco_sites):.1f}")
    
    # Physical properties
    pIs = [b["validation"]["properties"]["pI"] for b in binders]
    gravys = [b["validation"]["properties"]["GRAVY"] for b in binders]
    print(f"\nPhysicochemical Properties:")
    print(f"pI range: {min(pIs):.1f} - {max(pIs):.1f}")
    print(f"GRAVY range: {min(gravys):.3f} - {max(gravys):.3f}")
    
    # Sequence similarity groups
    print(f"\nSequence Similarity:")
    if results["similar_groups"]:
        print(f"Found {len(results['similar_groups'])} groups of similar sequences")
        for i, group in enumerate(results['similar_groups'], 1):
            print(f"Group {i}: sequences {', '.join(str(idx) for idx in group)}")
    else:
        print("No highly similar sequences found")

if __name__ == "__main__":
    input_file = "output/codette_antibody_designs_20250912_150658.json"
    output_file = f"output/validation_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    print(f"Running validation on {input_file}...")
    results = validate_binder_set(input_file, output_file)
    
    print(f"\nValidation complete. Results saved to {output_file}")
    analyze_validation_results(results)