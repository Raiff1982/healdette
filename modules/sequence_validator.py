"""
Sequence validation module for comparing generated antibodies against known therapeutic antibodies.
"""

import json
import os
from typing import Dict, List, Tuple
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class TherapeuticAntibodyValidator:
    def __init__(self, validation_data_path: str = None):
        """
        Initialize validator with therapeutic antibody dataset.
        
        Args:
            validation_data_path: Path to validation dataset JSON file
        """
        if validation_data_path is None:
            validation_data_path = os.path.join(
                os.path.dirname(__file__), 
                'data', 
                'therapeutic_antibodies.json'
            )
        
        with open(validation_data_path, 'r') as f:
            self.validation_data = json.load(f)
        
        # Extract validation metrics and data
        self.metrics = self.validation_data['validation_metrics']
        self.antibodies = self.validation_data['therapeutic_antibodies']
        
        # Default configuration
        self.default_config = {
            "sequence_properties": {
                "min_length": 100,
                "max_length": 500,
                "max_hydrophobicity": 0.5,
                "max_instability": 40.0,
                "min_antigenic_regions": 2
            },
            "structure_requirements": {
                "min_helix_content": 0.2,
                "min_sheet_content": 0.1,
                "max_disorder_regions": 3
            },
            "validation_thresholds": {
                "similarity_cutoff": 0.7,
                "quality_score_min": 0.5,
                "confidence_threshold": 0.8
            },
            "signal_peptide": {
                "enabled": True,
                "min_length": 15,
                "max_length": 30,
                "required": False,
                "strip": False,
                "confidence_threshold": 0.6,
                "n_region_basic_threshold": 0.3,
                "h_region_hydrophobic_threshold": 0.6
            }
        }
        
        # Initialize configuration
        if not hasattr(self, 'config'):
            self.config = {}
        
    def extract_cdrs(self, sequence: str) -> Dict[str, List[str]]:
        """
        Extract CDRs from a sequence based on framework region patterns.
        
        Args:
            sequence: Full antibody sequence
            
        Returns:
            Dictionary containing CDR sequences
        """
        # Framework region patterns (simplified version)
        fr_patterns = {
            'FR1': 'EVQLVESGGGLVQPGGSLRLSCAAS',
            'FR2': 'WVRQAPGKGLEWV',
            'FR3': 'RFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
            'FR4': 'WGQGTLVTVSS'
        }
        
        cdrs = {
            'cdr1': '',
            'cdr2': '',
            'cdr3': ''
        }
        
        try:
            # Find FR1
            fr1_end = sequence.find(fr_patterns['FR1']) + len(fr_patterns['FR1'])
            # Find FR2
            fr2_start = sequence.find(fr_patterns['FR2'])
            # Find FR3
            fr3_start = sequence.find(fr_patterns['FR3'])
            # Find FR4
            fr4_start = sequence.find(fr_patterns['FR4'])
            
            if all(x != -1 for x in [fr1_end, fr2_start, fr3_start, fr4_start]):
                cdrs['cdr1'] = sequence[fr1_end:fr2_start]
                cdrs['cdr2'] = sequence[fr2_start + len(fr_patterns['FR2']):fr3_start]
                cdrs['cdr3'] = sequence[fr3_start + len(fr_patterns['FR3']):fr4_start]
        except Exception:
            # If pattern matching fails, return empty CDRs
            pass
            
        return cdrs
        
    def calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence similarity using position-specific scoring.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Similarity score between 0 and 1
        """
        # Define groups of similar amino acids
        similar_groups = [
            set('ILVM'),           # Aliphatic
            set('FYW'),           # Aromatic
            set('KRH'),           # Basic
            set('DE'),            # Acidic
            set('STNQ'),          # Polar
            set('AG'),            # Small
            set('C'),             # Cysteine
            set('P')              # Proline
        ]
        
        def aa_similarity(aa1: str, aa2: str) -> float:
            """Calculate similarity score between two amino acids."""
            if aa1 == aa2:
                return 1.0
            for group in similar_groups:
                if aa1 in group and aa2 in group:
                    return 0.5
            return 0.0
        
        # Get the length of the shorter sequence
        min_length = min(len(seq1), len(seq2))
        
        # Calculate similarity scores for the overlapping region
        similarity_scores = [
            aa_similarity(aa1, aa2) 
            for aa1, aa2 in zip(seq1[:min_length], seq2[:min_length])
        ]
        
        # Add penalties for length difference
        length_difference = abs(len(seq1) - len(seq2))
        length_penalty = length_difference * 0.1
        
        # Calculate final score
        if not similarity_scores:
            return 0.0
            
        raw_score = sum(similarity_scores) / len(similarity_scores)
        final_score = max(0.0, min(1.0, raw_score - length_penalty))
        
        return final_score
        
    def find_similar_antibodies(self, sequence: Dict[str, str], threshold: float = 0.7) -> List[Dict]:
        """
        Find similar therapeutic antibodies in the validation set.
        
        Args:
            sequence: Dictionary with heavy_chain and light_chain sequences
            threshold: Minimum similarity score to consider
            
        Returns:
            List of similar antibodies with similarity scores
        """
        similar_antibodies = []
        
        for antibody in self.antibodies:
            heavy_similarity = self.calculate_sequence_similarity(
                sequence['heavy_chain'],
                antibody['sequence']['heavy_chain']
            )
            
            light_similarity = self.calculate_sequence_similarity(
                sequence['light_chain'],
                antibody['sequence']['light_chain']
            )
            
            # Average similarity score weighted slightly towards heavy chain
            overall_similarity = (0.6 * heavy_similarity + 0.4 * light_similarity)
            
            if overall_similarity >= threshold:
                similar_antibodies.append({
                    'name': antibody['name'],
                    'similarity_score': round(overall_similarity, 3),
                    'heavy_chain_similarity': round(heavy_similarity, 3),
                    'light_chain_similarity': round(light_similarity, 3),
                    'target': antibody['target'],
                    'antibody_type': antibody['antibody_type']
                })
                
        return sorted(similar_antibodies, key=lambda x: x['similarity_score'], reverse=True)
        
    def validate_sequence_metrics(self, sequence: Dict[str, str]) -> Dict:
        """
        Validate sequence metrics against known therapeutic antibodies.
        
        Args:
            sequence: Dictionary with heavy_chain and light_chain sequences
            
        Returns:
            Dictionary containing validation results
        """
        results = {
            'valid': True,
            'warnings': [],
            'metrics': {}
        }
        
        # Check sequence lengths
        for chain in ['heavy_chain', 'light_chain']:
            length = len(sequence[chain])
            valid_range = self.metrics['sequence_length'][chain]
            
            results['metrics'][f'{chain}_length'] = length
            
            if length < valid_range['min'] or length > valid_range['max']:
                results['valid'] = False
                results['warnings'].append(
                    f'{chain} length {length} outside valid range '
                    f'({valid_range["min"]}-{valid_range["max"]})'
                )
        
        # Extract and validate CDRs
        for chain in ['heavy_chain', 'light_chain']:
            cdrs = self.extract_cdrs(sequence[chain])
            results['metrics'][f'{chain}_cdrs'] = cdrs
            
            for cdr_name, cdr_seq in cdrs.items():
                if cdr_seq:  # Only validate if CDR was found
                    valid_range = self.metrics['cdr_length_ranges'][chain.split('_')[0]][cdr_name]
                    
                    if len(cdr_seq) < valid_range[0] or len(cdr_seq) > valid_range[1]:
                        results['valid'] = False
                        results['warnings'].append(
                            f'{chain} {cdr_name} length {len(cdr_seq)} outside valid range '
                            f'({valid_range[0]}-{valid_range[1]})'
                        )
                else:
                    results['valid'] = False
                    results['warnings'].append(f'Could not identify {cdr_name} in {chain}')
        
        # Find similar antibodies
        similar = self.find_similar_antibodies(sequence)
        if similar:
            results['similar_antibodies'] = similar
            
            # Warning if too similar to existing antibody
            if similar[0]['similarity_score'] > 0.9:
                results['warnings'].append(
                    f'Very high similarity ({similar[0]["similarity_score"]:.3f}) to '
                    f'existing antibody {similar[0]["name"]}'
                )
        
        return results

def validate_generated_sequences(sequences: List[Dict[str, str]], 
                              output_file: str = None) -> Dict:
    """
    Validate a set of generated sequences against therapeutic antibody dataset.
    
    Args:
        sequences: List of sequences to validate
        output_file: Optional path to save validation results
        
    Returns:
        Dictionary containing validation results
    """
    validator = TherapeuticAntibodyValidator()
    
    results = {
        'validated_sequences': [],
        'summary': {
            'total': len(sequences),
            'valid': 0,
            'warnings': 0,
            'similar_to_existing': 0
        }
    }
    
    for sequence in sequences:
        validation = validator.validate_sequence_metrics(sequence)
        
        results['validated_sequences'].append({
            'sequence': sequence,
            'validation': validation
        })
        
        # Update summary statistics
        if validation['valid']:
            results['summary']['valid'] += 1
        if validation['warnings']:
            results['summary']['warnings'] += 1
        if validation.get('similar_antibodies'):
            results['summary']['similar_to_existing'] += 1
    
    # Calculate percentages
    total = results['summary']['total']
    if total > 0:
        results['summary']['valid_percentage'] = round(
            100 * results['summary']['valid'] / total, 1
        )
        results['summary']['warning_percentage'] = round(
            100 * results['summary']['warnings'] / total, 1
        )
        results['summary']['similarity_percentage'] = round(
            100 * results['summary']['similar_to_existing'] / total, 1
        )
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
    
    return results