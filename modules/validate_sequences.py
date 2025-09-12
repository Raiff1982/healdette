"""
Comprehensive validation module for antibody sequences.
Performs computational checks for various sequence properties and potential issues.
"""

import math
import numpy as np
from collections import defaultdict
import re
import json
from typing import Dict, List, Tuple

class SequenceValidator:
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        
    def predict_disorder(self) -> float:
        """
        Simple disorder prediction based on amino acid propensities.
        Returns fraction of residues predicted to be disordered.
        """
        # Disorder-promoting residues (based on literature)
        disorder_prone = set('RKEPNDQSG')
        disorder_count = sum(1 for aa in self.sequence if aa in disorder_prone)
        return disorder_count / len(self.sequence)
    
    def check_signal_peptide(self) -> Dict:
        """
        Basic signal peptide detection heuristics.
        Returns likelihood of being a signal peptide.
        """
        if len(self.sequence) < 15:
            return {"has_signal": False, "confidence": 1.0}
            
        n_region = self.sequence[:5]
        h_region = self.sequence[5:12]
        c_region = self.sequence[12:15]
        
        # Check for typical signal peptide features
        n_region_positive = sum(aa in 'KR' for aa in n_region) >= 1
        h_region_hydrophobic = sum(aa in 'AILMFWV' for aa in h_region) >= 4
        c_region_pattern = bool(re.search('[AGST].[AGST]', c_region))
        
        confidence = (n_region_positive + h_region_hydrophobic + c_region_pattern) / 3
        return {
            "has_signal": confidence > 0.6,
            "confidence": confidence
        }
    
    def analyze_cysteines(self) -> Dict:
        """
        Analyze cysteine content and potential disulfide bonds.
        """
        cys_positions = [i for i, aa in enumerate(self.sequence) if aa == 'C']
        n_cys = len(cys_positions)
        
        return {
            "count": n_cys,
            "paired": n_cys % 2 == 0,
            "positions": cys_positions,
            "spacing": [cys_positions[i+1] - cys_positions[i] for i in range(len(cys_positions)-1)] if n_cys > 1 else []
        }
    
    def find_glycosylation_sites(self) -> List[Dict]:
        """
        Identify potential N-glycosylation sites (N-X-S/T).
        """
        pattern = re.compile('N[^P][ST]')
        sites = []
        
        for match in pattern.finditer(self.sequence):
            sites.append({
                "position": match.start(),
                "motif": self.sequence[match.start():match.start()+3]
            })
        
        return sites
    
    def calculate_properties(self) -> Dict:
        """
        Calculate various physicochemical properties.
        """
        # Kyte & Doolittle hydropathy values
        hydropathy = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        # Calculate GRAVY (Grand Average of Hydropathy)
        gravy = sum(hydropathy[aa] for aa in self.sequence) / len(self.sequence)
        
        # Calculate molecular weight
        weights = {
            'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
            'Q': 146.2, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
            'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
            'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
        }
        mw = sum(weights[aa] for aa in self.sequence)
        
        # pKa values for amino acids
        pka_values = {
            'K': 10.0,  # Lysine
            'R': 12.0,  # Arginine
            'H': 6.08,  # Histidine
            'D': 3.65,  # Aspartic acid
            'E': 4.25,  # Glutamic acid
            'C': 8.18,  # Cysteine
            'Y': 10.46, # Tyrosine
            'N_term': 8.0,  # N-terminus
            'C_term': 3.1   # C-terminus
        }
        
        # Calculate pI using iterative method
        def charge_at_ph(ph):
            charge = 0
            # N-terminus
            charge += 1 / (1 + 10**(ph - pka_values['N_term']))
            # C-terminus
            charge -= 1 / (1 + 10**(pka_values['C_term'] - ph))
            
            for aa in self.sequence:
                if aa in 'KRH':  # Basic residues
                    charge += 1 / (1 + 10**(ph - pka_values[aa]))
                elif aa in 'DE':  # Acidic residues
                    charge -= 1 / (1 + 10**(pka_values[aa] - ph))
                elif aa == 'C':  # Cysteine
                    charge -= 1 / (1 + 10**(pka_values[aa] - ph))
                elif aa == 'Y':  # Tyrosine
                    charge -= 1 / (1 + 10**(pka_values[aa] - ph))
            return charge
        
        # Binary search for pI
        ph_min, ph_max = 0, 14
        while ph_max - ph_min > 0.01:
            ph_mid = (ph_min + ph_max) / 2
            charge = charge_at_ph(ph_mid)
            if charge > 0:
                ph_min = ph_mid
            else:
                ph_max = ph_mid
        
        return {
            "pI": round((ph_min + ph_max) / 2, 2),
            "GRAVY": gravy,
            "molecular_weight": mw,
            "aromaticity": sum(aa in 'FWY' for aa in self.sequence) / len(self.sequence),
            "instability_index": None  # Would need complex calculation
        }
    
    @staticmethod
    def calculate_similarity(seq1: str, seq2: str) -> float:
        """
        Calculate sequence similarity between two sequences.
        """
        if len(seq1) != len(seq2):
            return 0.0
        matches = sum(a == b for a, b in zip(seq1, seq2))
        return matches / len(seq1)

def validate_binder(sequence: str) -> Dict:
    """
    Perform comprehensive validation of a single binder sequence.
    """
    validator = SequenceValidator(sequence)
    
    return {
        "length": len(sequence),
        "disorder": validator.predict_disorder(),
        "signal_peptide": validator.check_signal_peptide(),
        "cysteines": validator.analyze_cysteines(),
        "glycosylation": validator.find_glycosylation_sites(),
        "properties": validator.calculate_properties()
    }

def validate_binder_set(json_file: str, output_file: str = None):
    """
    Validate a set of binders from a JSON file and optionally save results.
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    results = []
    for binder in data['personalized_binders']:
        validation = validate_binder(binder['sequence'])
        results.append({
            **binder,
            "validation": validation
        })
    
    # Group similar sequences
    similar_groups = []
    used = set()
    
    for i, binder1 in enumerate(results):
        if i in used:
            continue
        
        group = [i]
        for j, binder2 in enumerate(results[i+1:], i+1):
            if j not in used and SequenceValidator.calculate_similarity(
                binder1['sequence'], binder2['sequence']) > 0.9:
                group.append(j)
                used.add(j)
        
        if len(group) > 1:
            similar_groups.append(group)
    
    output = {
        "validated_binders": results,
        "similar_groups": similar_groups
    }
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=4)
    
    return output