from modules.generate_binders import AntibodyGenerator

def main():
    # Initialize generator
    generator = AntibodyGenerator()
    
    # Test context with a real protein motif
    test_context = {
        "cleaned_sequence": "QVQLVQSGAEVKKPGASVKVSCKASGYTFT",  # Common VH framework motif
    }
    
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
        print("\nFull Sequence:")
        print(binder['sequence'])
        print("-" * 50)
    
    print(f"\nGeneration Statistics:")
    print(f"Success Rate: {results['stats']['success_rate']:.2%}")
    print(f"Total Attempts: {results['stats']['attempts']}")

if __name__ == "__main__":
    main()