"""Sequence quality validation for antibody generation."""
from typing import Dict
import re

class SequenceValidator:
    """Validates antibody sequences for quality and composition."""
    
    def __init__(self):
        """Initialize sequence validator with quality parameters."""
        self.params = {
            'max_homopolymer_length': 4,
            'min_cdr_length': 5,
            'max_cdr_length': 20,
            'max_aa_frequency': 0.35,  # 35% max for any amino acid
            'min_unique_aa': 5,  # At least 5 different amino acids
            'hydrophobic_range': (0.15, 0.65),  # 15-65% hydrophobic content
            'hydrophobic_aas': set('AILMFWYV')
        }
    
    def validate_cdr(self, sequence: str) -> bool:
        """Validate CDR sequence meets all quality criteria."""
        # Length check
        if not (self.params['min_cdr_length'] <= len(sequence) <= self.params['max_cdr_length']):
            return False
        
        # Homopolymer check using regex
        for aa in set(sequence):
            pattern = f"{aa}{{{self.params['max_homopolymer_length']},}}"
            if re.search(pattern, sequence):
                return False
        
        # Amino acid composition checks
        aa_counts = {}
        for aa in sequence:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        # Check amino acid diversity
        if len(aa_counts) < self.params['min_unique_aa']:
            return False
        
        # Check maximum frequency
        sequence_length = len(sequence)
        for count in aa_counts.values():
            if count / sequence_length > self.params['max_aa_frequency']:
                return False
        
        # Hydrophobicity check
        hydrophobic_count = sum(1 for aa in sequence if aa in self.params['hydrophobic_aas'])
        hydrophobic_fraction = hydrophobic_count / sequence_length
        min_hydrophobic, max_hydrophobic = self.params['hydrophobic_range']
        
        if not (min_hydrophobic <= hydrophobic_fraction <= max_hydrophobic):
            return False
        
        return True
    
    def analyze_sequence(self, sequence: str) -> Dict:
        """Analyze sequence properties for detailed validation."""
        aa_counts = {}
        for aa in sequence:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        sequence_length = len(sequence)
        hydrophobic_count = sum(1 for aa in sequence if aa in self.params['hydrophobic_aas'])
        
        # Find longest homopolymer
        max_homopolymer = 1
        current_run = 1
        prev_aa = sequence[0]
        
        for aa in sequence[1:]:
            if aa == prev_aa:
                current_run += 1
                max_homopolymer = max(max_homopolymer, current_run)
            else:
                current_run = 1
            prev_aa = aa
        
        return {
            'length': sequence_length,
            'unique_aa_count': len(aa_counts),
            'max_aa_frequency': max(count / sequence_length for count in aa_counts.values()),
            'most_common_aa': max(aa_counts.items(), key=lambda x: x[1])[0],
            'hydrophobic_fraction': hydrophobic_count / sequence_length,
            'max_homopolymer': max_homopolymer,
            'passes_validation': self.validate_cdr(sequence)
        }