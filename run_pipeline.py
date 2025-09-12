"""
Reproducibility harness for Healdette pipeline.
Runs the full pipeline with validation and generates all artifacts.
"""

import argparse
import json
import sys
import os
import hashlib
from datetime import datetime
import numpy as np
import torch
import pandas as pd

from modules.validate_sequences import validate_binder_set

def set_random_seeds(seed=42):
    """Set random seeds for reproducibility."""
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)

def calculate_sha256(filepath):
    """Calculate SHA256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def validate_criteria(results, criteria):
    """Validate results against pre-registered criteria."""
    failures = []
    for binder in results['validated_binders']:
        validation = binder['validation']
        
        # Check disorder
        if validation['disorder'] > criteria['disorder_threshold']:
            failures.append(f"Sequence {binder['sequence'][:20]}... has high disorder: {validation['disorder']:.3f}")
            
        # Check signal peptide
        if criteria['signal_peptide'] == 'disallow' and validation['signal_peptide']['has_signal']:
            failures.append(f"Sequence {binder['sequence'][:20]}... has signal peptide")
            
        # Check cysteine pairs
        if criteria['cys_pairs'] == 'required' and not validation['cysteines']['patterns']['paired']:
            failures.append(f"Sequence {binder['sequence'][:20]}... lacks paired cysteines")
            
        # Check GRAVY
        gravy = validation['properties']['GRAVY']
        if not (criteria['gravy_range'][0] <= gravy <= criteria['gravy_range'][1]):
            failures.append(f"Sequence {binder['sequence'][:20]}... has GRAVY {gravy:.3f} outside range")
    
    return failures

def generate_triage_table(results):
    """Generate triage table with key metrics."""
    rows = []
    for binder in results['validated_binders']:
        rows.append({
            'sequence_length': len(binder['sequence']),
            'personalization_score': binder['personalization_score'],
            'disorder': binder['validation']['disorder'],
            'cys_pairs': binder['validation']['cysteines']['count'] // 2,
            'glyco_sites': len(binder['validation']['glycosylation']),
            'gravy': binder['validation']['properties']['GRAVY'],
            'pI': binder['validation']['properties']['pI']
        })
    
    df = pd.DataFrame(rows)
    return df.sort_values('personalization_score', ascending=False)

def main(args):
    # Load configuration
    with open('run_manifest.json', 'r') as f:
        manifest = json.load(f)
    
    # Set deterministic mode if requested
    if args.deterministic:
        set_random_seeds()
    
    # Run validation
    results = validate_binder_set(args.input_json)
    
    # Generate triage table
    triage_table = generate_triage_table(results)
    triage_table.to_csv('output/triage_table.csv')
    
    # Validate against criteria
    failures = validate_criteria(results, manifest['validation_criteria'])
    
    # Calculate checksums
    checksums = {}
    for filepath in [args.input_json, 'output/triage_table.csv', 'output/sequence_analysis.png']:
        checksums[os.path.basename(filepath)] = calculate_sha256(filepath)
    
    with open('checksums.sha256', 'w') as f:
        for filename, checksum in checksums.items():
            f.write(f"{checksum}  {filename}\n")
    
    # Exit with error if validation failed
    if failures:
        print("\nValidation failures:")
        for failure in failures:
            print(f"- {failure}")
        sys.exit(1)
    
    print("\nValidation successful!")
    print(f"Results saved to {args.output_dir}")
    sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Healdette pipeline with validation')
    parser.add_argument('--input-json', default='output/codette_antibody_designs_20250912_150658.json',
                      help='Input JSON file with antibody designs')
    parser.add_argument('--output-dir', default='output',
                      help='Output directory for results')
    parser.add_argument('--deterministic', action='store_true',
                      help='Run in deterministic mode with fixed seeds')
    
    args = parser.parse_args()
    main(args)