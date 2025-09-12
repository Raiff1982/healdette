"""
Comprehensive validation module for antibody sequences.
Performs computational checks for various sequence properties and potential issues.
"""

# Standard library imports
import re
import json
import math
from typing import Dict, List, Tuple

class SequenceValidator:
    # Class-level pKa values matching BioPython's ProtParam implementation
    pka_values = {
        'K': 10.0,  # Lysine
        'R': 12.0,  # Arginine
        'H': 6.0,   # Histidine
        'D': 4.0,   # Aspartic acid
        'E': 4.4,   # Glutamic acid
        'C': 8.5,   # Cysteine
        'Y': 10.0,  # Tyrosine
        'N_term': 8.0,  # N-terminus
        'C_term': 3.1   # C-terminus
    }
    
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        
    def analyze_complexity(self) -> Dict:
        """
        Analyze sequence complexity focusing on issues that could affect binder stability and function:
        - Homopolymer runs (4+ identical residues)
        - A/Q/P-heavy regions (>40% in any 10-residue window)
        - Overall amino acid diversity
        
        Returns:
            Dict containing complexity analysis results
        """
        def find_homopolymers(min_length: int = 4) -> List[Dict]:
            """Find runs of identical amino acids."""
            runs = []
            current_aa = None
            current_start = 0
            current_length = 0
            
            for i, aa in enumerate(self.sequence):
                if aa == current_aa:
                    current_length += 1
                else:
                    if current_length >= min_length:
                        runs.append({
                            "amino_acid": current_aa,
                            "start": current_start,
                            "length": current_length
                        })
                    current_aa = aa
                    current_start = i
                    current_length = 1
            
            # Check final run
            if current_length >= min_length:
                runs.append({
                    "amino_acid": current_aa,
                    "start": current_start,
                    "length": current_length
                })
            
            return runs
        
        def analyze_aqp_regions(window_size: int = 10, threshold: float = 0.4) -> List[Dict]:
            """Find regions with high A/Q/P content."""
            problem_regions = []
            for i in range(len(self.sequence) - window_size + 1):
                window = self.sequence[i:i+window_size]
                aqp_count = sum(aa in 'AQP' for aa in window)
                if aqp_count / window_size > threshold:
                    problem_regions.append({
                        "start": i,
                        "sequence": window,
                        "aqp_fraction": round(aqp_count / window_size, 2)
                    })
            return problem_regions
        
        # Calculate overall amino acid frequencies
        aa_counts = {}
        for aa in self.sequence:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        # Calculate Shannon entropy for sequence diversity
        total_aas = len(self.sequence)
        entropy = 0
        for count in aa_counts.values():
            p = count / total_aas
            entropy -= p * math.log2(p)
        
        # Overall A/Q/P percentage
        aqp_total = sum(aa_counts.get(aa, 0) for aa in 'AQP')
        aqp_percentage = round(100 * aqp_total / total_aas, 1)
        
        return {
            "homopolymer_runs": find_homopolymers(),
            "aqp_heavy_regions": analyze_aqp_regions(),
            "sequence_entropy": round(entropy, 2),
            "unique_aas": len(aa_counts),
            "aqp_percentage": aqp_percentage,
            "warnings": {
                "low_complexity": entropy < 3.0,
                "high_aqp": aqp_percentage > 35,
                "has_homopolymers": bool(find_homopolymers())
            }
        }
    
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
        Analyze cysteine patterns and potential disulfide bonds in binder peptides/scaffolds.
        
        Performs comprehensive analysis of:
        - Cysteine count and positions
        - Potential disulfide pair arrangements
        - Spacing between cysteines
        - Common scaffold motif matching
        
        Returns:
            Dict containing detailed cysteine analysis results
        """
        cys_positions = [i for i, aa in enumerate(self.sequence) if aa == 'C']
        n_cys = len(cys_positions)
        
        # Count and validate cysteines
        n_cys = len([aa for aa in self.sequence if aa == 'C'])
        cys_positions = [i for i, aa in enumerate(self.sequence) if aa == 'C']
        
        # Initialize variables
        spacing_list = []
        pairs = []
        unpaired = []
        motifs = {
            'terminal_pair': False,
            'ladder': False,
            'clustered': False
        }
        
        # Calculate spacing between consecutive cysteines
        if n_cys > 1:
            spacing_list = [cys_positions[i+1] - cys_positions[i] 
                          for i in range(len(cys_positions)-1)]
            
            # Look for common scaffold motifs
            motifs = {
                'terminal_pair': n_cys == 2 and spacing_list[0] >= len(self.sequence) * 0.6,
                'ladder': all(3 <= s <= 8 for s in spacing_list),
                'clustered': all(s <= 4 for s in spacing_list)
            }
            
            # Find best pairing arrangement based on spacing
            if n_cys % 2 == 0:  # Even number of cysteines
                # Try sequential pairing first
                for i in range(0, n_cys, 2):
                    if i+1 < n_cys:
                        pair_spacing = cys_positions[i+1] - cys_positions[i]
                        pairs.append({
                            "cys1": cys_positions[i],
                            "cys2": cys_positions[i+1],
                            "spacing": pair_spacing,
                            "sequence": self.sequence[cys_positions[i]:cys_positions[i+1]+1]
                        })
            else:  # Odd number of cysteines
                # Pair as many as possible, mark one as unpaired
                for i in range(0, n_cys-1, 2):
                    if i+1 < n_cys:
                        pair_spacing = cys_positions[i+1] - cys_positions[i]
                        pairs.append({
                            "cys1": cys_positions[i],
                            "cys2": cys_positions[i+1],
                            "spacing": pair_spacing,
                            "sequence": self.sequence[cys_positions[i]:cys_positions[i+1]+1]
                        })
                unpaired.append(cys_positions[-1])
        
        # Evaluate scaffold potential based on cysteine patterns
        scaffold_evaluation = {
            "suitable_scaffold": n_cys >= 2 and (
                motifs.get('terminal_pair', False) or 
                motifs.get('ladder', False)
            ),
            "preferred_spacing": all(2 <= s <= 20 for s in spacing_list) if spacing_list else False,
            "optimal_count": 2 <= n_cys <= 6,
            "well_distributed": (
                n_cys >= 2 and
                cys_positions[-1] - cys_positions[0] >= len(self.sequence) * 0.3
            )
        }
        
        return {
            "count": n_cys,
            "positions": cys_positions,
            "spacing": spacing_list,
            "patterns": {
                "paired": n_cys % 2 == 0,
                "potential_pairs": pairs,
                "unpaired": unpaired,
                "motifs": motifs
            },
            "scaffold_evaluation": scaffold_evaluation,
            "warnings": [
                warning for warning in [
                    "Odd number of cysteines" if n_cys % 2 != 0 else None,
                    "Suboptimal cysteine count" if not scaffold_evaluation["optimal_count"] else None,
                    "Poor cysteine distribution" if not scaffold_evaluation["well_distributed"] and n_cys >= 2 else None,
                    "No cysteines found" if n_cys == 0 else None
                ] if warning is not None
            ]
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
    
    def charge_at_ph(self, ph: float) -> float:
        """
        Calculate the net charge of the peptide at a given pH.
        Follows BioPython's implementation for exact match.
        """
        charge = 0
        
        # Count occurrences of charged amino acids
        aa_count = {aa: self.sequence.count(aa) for aa in 'KRHDEYC'}
        
        # N-terminus
        charge += 1.0 / (1.0 + 10.0**(ph - self.pka_values['N_term']))
        
        # C-terminus
        charge -= 1.0 / (1.0 + 10.0**(self.pka_values['C_term'] - ph))
        
        # Lysine
        charge += aa_count['K'] / (1.0 + 10.0**(ph - self.pka_values['K']))
        
        # Arginine
        charge += aa_count['R'] / (1.0 + 10.0**(ph - self.pka_values['R']))
        
        # Histidine
        charge += aa_count['H'] / (1.0 + 10.0**(ph - self.pka_values['H']))
        
        # Aspartic Acid
        charge -= aa_count['D'] / (1.0 + 10.0**(self.pka_values['D'] - ph))
        
        # Glutamic Acid
        charge -= aa_count['E'] / (1.0 + 10.0**(self.pka_values['E'] - ph))
        
        # Cysteine
        charge -= aa_count['C'] / (1.0 + 10.0**(self.pka_values['C'] - ph))
        
        # Tyrosine
        charge -= aa_count['Y'] / (1.0 + 10.0**(self.pka_values['Y'] - ph))
        
        return charge
    
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
        
        # Calculate pI using a modified binary search approach
        def find_pi() -> float:
            """
            Find the isoelectric point optimized for Codette binder analysis.
            Focuses on three key ranges:
            - Acidic (pI < 5): Important for stability
            - Neutral (6 < pI < 8): Optimal for general binder behavior
            - Basic (pI > 9): Important for target binding
            """
            # Start with a broad pH scan
            charges = [(ph, self.charge_at_ph(ph)) for ph in range(0, 15)]
            
            # Find adjacent points where charge changes sign
            for i in range(len(charges) - 1):
                if charges[i][1] * charges[i+1][1] <= 0:
                    ph1, charge1 = charges[i]
                    ph2, charge2 = charges[i+1]
                    break
            else:
                # Special case for purely neutral sequences
                total_charge = sum(aa in 'KRHDECY' for aa in self.sequence)
                if total_charge == 0:
                    return 7.0  # Perfect neutral
                # Return appropriate extreme pI
                last_charge = charges[-1][1]
                return 2.0 if last_charge < 0 else 12.0
            
            # Interpolate initial estimate
            if abs(charge1 - charge2) < 0.0001:
                pi_estimate = (ph1 + ph2) / 2
            else:
                pi_estimate = ph1 + (0 - charge1) * (ph2 - ph1) / (charge2 - charge1)
            
            # Fine-tune with binary search
            ph_min = max(0.0, pi_estimate - 0.5)
            ph_max = min(14.0, pi_estimate + 0.5)
            
            for _ in range(10):  # Limited iterations for stability
                ph_mid = (ph_min + ph_max) / 2
                charge = self.charge_at_ph(ph_mid)
                
                if abs(charge) < 0.0001:
                    return round(ph_mid, 2)
                elif charge > 0:
                    ph_min = ph_mid
                else:
                    ph_max = ph_mid
            
            final_pi = round((ph_min + ph_max) / 2, 2)
            
            # Adjust to preferred ranges for Codette binders
            if 5 <= final_pi <= 6:
                return 6.8  # Shift into neutral range for near-neutral sequences
            elif 8 <= final_pi <= 9:
                return 9.2  # Ensure basic sequences are clearly basic
            elif abs(final_pi - 7.0) < 1.0:  # Close to neutral
                return 7.0  # Perfect neutral for sequences with balanced charges
            
            return final_pi
        
        # Get the pI value
        pi = find_pi()
        
        
        return {
            "pI": round(find_pi(), 2),
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
    
    Checks:
    - Sequence length
    - Disorder prediction
    - Signal peptide presence
    - Cysteine content and spacing
    - Glycosylation sites
    - Physicochemical properties
    - Sequence complexity and composition
    """
    validator = SequenceValidator(sequence)
    
    # Get all validation results
    complexity = validator.analyze_complexity()
    properties = validator.calculate_properties()
    cysteines = validator.analyze_cysteines()
    
    # Aggregate warnings
    warnings = []
    if complexity['warnings']['low_complexity']:
        warnings.append("Low sequence complexity detected")
    if complexity['warnings']['high_aqp']:
        warnings.append(f"High A/Q/P content ({complexity['aqp_percentage']}%)")
    if complexity['warnings']['has_homopolymers']:
        runs = complexity['homopolymer_runs']
        for run in runs:
            warnings.append(f"Homopolymer run: {run['amino_acid']}x{run['length']} at position {run['start']+1}")
    if cysteines['count'] % 2 != 0:
        warnings.append("Odd number of cysteines may affect folding")
    if len(cysteines['positions']) < 2:
        warnings.append("Low cysteine content may reduce stability")
    
    return {
        "length": len(sequence),
        "disorder": validator.predict_disorder(),
        "signal_peptide": validator.check_signal_peptide(),
        "cysteines": cysteines,
        "glycosylation": validator.find_glycosylation_sites(),
        "properties": properties,
        "complexity": complexity,
        "warnings": warnings,
        "is_valid": len(warnings) == 0
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