"""
Debug script to analyze pI calculation differences.
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from modules.validate_sequences import SequenceValidator

def debug_sequence(sequence: str):
    print(f"\nAnalyzing sequence: {sequence}")
    
    # Our implementation
    validator = SequenceValidator(sequence)
    our_pi = validator.calculate_properties()['pI']
    
    # BioPython's implementation
    biopython = ProteinAnalysis(sequence)
    bio_pi = round(biopython.isoelectric_point(), 2)
    
    print(f"Our pI: {our_pi}")
    print(f"BioPython pI: {bio_pi}")
    
    # Print charge values at various pH levels
    print("\npH\tOur Charge\tBioPython Charge")
    print("-" * 40)
    
    for ph in [2.0, 3.0, 4.0, 4.05, 4.1, 5.0, 6.0, 7.0]:
        our_charge = validator.charge_at_ph(ph)
        print(f"{ph:.2f}\t{our_charge:.3f}\t\t{our_charge:.3f}")

if __name__ == '__main__':
    sequence = 'DDDEEEDDDEEE'  # The acidic sequence that's causing issues
    debug_sequence(sequence)