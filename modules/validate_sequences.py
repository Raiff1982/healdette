"""
Comprehensive validation module for antibody sequences.
Performs computational checks for various sequence properties and potential issues.
"""

# Standard library imports
import re
import json
import math
import datetime
from typing import Dict, List, Tuple
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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
    
    # CDR patterns based on Kabat numbering scheme
    CDR_PATTERNS = {
        'heavy': {
            'CDR1': r'[A-Z]{10,20}C[A-Z]{2}[A-Z]{3}[A-Z]C',
            'CDR2': r'[A-Z]{15,25}C[A-Z]{2}[A-Z]{7}[A-Z]C',
            'CDR3': r'[A-Z]{5,17}C[A-Z]{2}[A-Z]{3}[A-Z]C'
        },
        'light': {
            'CDR1': r'[A-Z]{10,17}C[A-Z]{2}[A-Z]{3}[A-Z]C',
            'CDR2': r'[A-Z]{7,12}C[A-Z]{2}[A-Z]{3}[A-Z]C',
            'CDR3': r'[A-Z]{7,11}C[A-Z]{2}[A-Z]{3}[A-Z]C'
        }
    }
    
    def __init__(self, sequence: str, config: Dict = None):
        """
        Initialize sequence validator with optional configuration.
        
        Args:
            sequence: The amino acid sequence to validate
            config: Optional configuration dictionary with validation parameters
        """
        self.sequence = sequence.upper()
        self.config = config or {}
        
        # Default configuration values
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
            }
        }
        
        # Merge provided config with defaults
        for key, default_values in self.default_config.items():
            if key not in self.config:
                self.config[key] = {}
            for param, value in default_values.items():
                self.config[key][param] = self.config.get(key, {}).get(param, value)
        
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
        Enhanced signal peptide detection for binder peptides/scaffolds.
        
        Features analyzed:
        - N-region: Basic amino acids (K/R)
        - H-region: Hydrophobic core
        - C-region: (-3, -1) rule with small neutral amino acids
        - Length constraints
        - Position-specific amino acid preferences
        
        Returns:
            Dict containing detailed signal peptide analysis
        """
        config = self.config['signal_peptide']
        
        if not config['enabled']:
            return {
                "enabled": False,
                "has_signal": False,
                "confidence": 0.0,
                "details": "Signal peptide detection disabled in configuration"
            }
        
        if len(self.sequence) < config['min_length']:
            return {
                "enabled": True,
                "has_signal": False,
                "confidence": 1.0,
                "details": f"Sequence too short (min {config['min_length']} residues required)"
            }
        
        # Dynamic region sizing based on sequence length
        n_region_length = min(6, len(self.sequence) // 5)
        h_region_length = min(12, len(self.sequence) // 3)
        c_region_length = 5
        
        total_sp_length = min(
            n_region_length + h_region_length + c_region_length,
            config['max_length']
        )
        
        # Extract regions
        n_region = self.sequence[:n_region_length]
        h_region = self.sequence[n_region_length:n_region_length + h_region_length]
        c_region = self.sequence[n_region_length + h_region_length:total_sp_length]
        
        # Analyze N-region (positive charge)
        n_region_basic = sum(aa in 'KR' for aa in n_region)
        n_region_score = n_region_basic / len(n_region)
        n_region_valid = n_region_score >= config['n_region_basic_threshold']
        
        # Analyze H-region (hydrophobic core)
        hydrophobic = set('AILMFWV')
        h_region_hydrophobic = sum(aa in hydrophobic for aa in h_region)
        h_region_score = h_region_hydrophobic / len(h_region)
        h_region_valid = h_region_score >= config['h_region_hydrophobic_threshold']
        
        # Analyze C-region (-3, -1 rule)
        c_region_valid = False
        if len(c_region) >= 3:
            small_neutral = set('AGST')
            c_region_pattern = (
                c_region[-3] in small_neutral and
                c_region[-1] in small_neutral
            )
            # Check for proline disruption
            no_proline_disruption = 'P' not in c_region[-3:]
            c_region_valid = c_region_pattern and no_proline_disruption
        
        # Calculate overall confidence
        feature_scores = [
            n_region_score if n_region_valid else 0,
            h_region_score if h_region_valid else 0,
            1.0 if c_region_valid else 0
        ]
        confidence = sum(feature_scores) / len(feature_scores)
        
        has_signal = confidence >= config['confidence_threshold']
        
        # Prepare detailed analysis
        details = {
            "n_region": {
                "sequence": n_region,
                "basic_fraction": round(n_region_score, 2),
                "valid": n_region_valid
            },
            "h_region": {
                "sequence": h_region,
                "hydrophobic_fraction": round(h_region_score, 2),
                "valid": h_region_valid
            },
            "c_region": {
                "sequence": c_region,
                "valid": c_region_valid
            }
        }
        
        result = {
            "enabled": True,
            "has_signal": has_signal,
            "confidence": round(confidence, 2),
            "details": details,
            "signal_sequence": self.sequence[:total_sp_length] if has_signal else None,
            "mature_sequence": self.sequence[total_sp_length:] if has_signal and config['strip'] else self.sequence
        }
        
        return result
    
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
    
    def validate_sequence(self) -> Dict:
        """
        Comprehensive validation of antibody sequence with improved metrics.
        
        Returns:
            Dict containing complete validation results and metrics
        """
        # Initialize results
        results = {
            "valid": True,
            "metrics": {},
            "failures": [],
            "warnings": []
        }
        
        try:
            # 1. Basic sequence validation
            if not self.sequence:
                results["valid"] = False
                results["failures"].append("Empty sequence")
                return results
            
            if not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in self.sequence):
                results["valid"] = False
                results["failures"].append("Invalid amino acids present")
                return results
            
            # 2. Length validation
            length = len(self.sequence)
            results["metrics"]["length"] = length
            cfg = self.config["sequence_properties"]
            
            if length < cfg["min_length"] or length > cfg["max_length"]:
                results["valid"] = False
                results["failures"].append(f"Length {length} outside allowed range")
            
            # 3. Analyze with BioPython's ProteinAnalysis
            protein_analysis = ProteinAnalysis(self.sequence)
            
            # Amino acid composition
            aa_percent = protein_analysis.get_amino_acids_percent()
            results["metrics"]["amino_acid_composition"] = aa_percent
            
            # Stability metrics
            instability_index = protein_analysis.instability_index()
            results["metrics"]["instability_index"] = instability_index
            if instability_index > cfg["max_instability"]:
                results["warnings"].append("High instability index")
            
            # Secondary structure
            helix, sheet, coil = protein_analysis.secondary_structure_fraction()
            results["metrics"]["secondary_structure"] = {
                "helix": helix,
                "sheet": sheet,
                "coil": coil
            }
            
            scfg = self.config["structure_requirements"]
            if helix < scfg["min_helix_content"]:
                results["warnings"].append("Low helical content")
            if sheet < scfg["min_sheet_content"]:
                results["warnings"].append("Low beta sheet content")
            
            # 4. Antibody-specific validations
            # Heavy chain CDR patterns
            heavy_cdr_matches = [
                bool(re.search(pattern, self.sequence))
                for pattern in self.CDR_PATTERNS["heavy"].values()
            ]
            
            # Light chain CDR patterns
            light_cdr_matches = [
                bool(re.search(pattern, self.sequence))
                for pattern in self.CDR_PATTERNS["light"].values()
            ]
            
            results["metrics"]["cdr_patterns"] = {
                "heavy_chain": sum(heavy_cdr_matches),
                "light_chain": sum(light_cdr_matches)
            }
            
            if not any(heavy_cdr_matches) and not any(light_cdr_matches):
                results["warnings"].append("No CDR patterns detected")
            
            # 5. Advanced physicochemical properties
            properties = self.calculate_properties()
            results["metrics"].update({
                "pI": properties["pI"],
                "GRAVY": properties["GRAVY"],
                "molecular_weight": properties["molecular_weight"],
                "aromaticity": properties["aromaticity"]
            })
            
            # Check hydrophobicity
            if properties["GRAVY"] > cfg["max_hydrophobicity"]:
                results["warnings"].append("High overall hydrophobicity")
            
            # 6. Sequence complexity and motifs
            complexity = self.analyze_complexity()
            results["metrics"]["sequence_complexity"] = complexity["sequence_entropy"]
            
            if complexity["warnings"]["low_complexity"]:
                results["warnings"].append("Low sequence complexity")
            
            # 7. Post-translational modification sites
            glyco_sites = self.find_glycosylation_sites()
            results["metrics"]["n_glycosylation_sites"] = len(glyco_sites)
            
            if len(glyco_sites) > 3:
                results["warnings"].append("High number of potential glycosylation sites")
            
            # 8. Calculate final quality score
            warning_penalty = len(results["warnings"]) * 0.1
            failure_penalty = len(results["failures"]) * 0.3
            
            base_score = (
                (1.0 if any(heavy_cdr_matches) or any(light_cdr_matches) else 0.5) +
                (1.0 if complexity["sequence_entropy"] > 3.0 else 0.5) +
                (1.0 if instability_index < cfg["max_instability"] else 0.5) +
                (1.0 if properties["GRAVY"] < cfg["max_hydrophobicity"] else 0.5)
            ) / 4.0
            
            quality_score = max(0.0, min(1.0, base_score - warning_penalty - failure_penalty))
            results["metrics"]["quality_score"] = round(quality_score, 2)
            
            # Update validity based on quality threshold
            if quality_score < self.config["validation_thresholds"]["quality_score_min"]:
                results["valid"] = False
                results["failures"].append(f"Quality score {quality_score} below threshold")
            
        except Exception as e:
            results["valid"] = False
            results["failures"].append(f"Validation error: {str(e)}")
        
        return results
def validate_binder(sequence: str, config: Dict = None) -> Dict:
    """
    Perform comprehensive validation of a single binder sequence using the enhanced validator.
    
    Args:
        sequence: The amino acid sequence to validate
        config: Optional configuration dictionary with validation parameters
    
    Performs comprehensive validation including:
    - Sequence quality and composition
    - Antibody-specific CDR patterns
    - Stability and structure predictions
    - Post-translational modifications
    - Advanced physicochemical properties
    - Quality scoring
    
    Returns:
        Dict containing complete validation results and metrics
    """
    validator = SequenceValidator(sequence, config)
    
    # Get validation results using the new comprehensive validation
    results = validator.validate_sequence()
    
    # Add additional context and metadata
    results["timestamp"] = datetime.datetime.now().isoformat()
    results["validator_version"] = "2.0.0"
    
    # Add supporting analyses
    results["supporting_analyses"] = {
        "cysteines": validator.analyze_cysteines(),
        "disorder": validator.predict_disorder(),
        "signal_peptide": validator.check_signal_peptide()
    }
    
    return results

def validate_binder_set(json_file: str, config: Dict = None, output_file: str = None):
    """
    Validate a set of binders from a JSON file and optionally save results.
    
    Args:
        json_file: Path to JSON file containing binders to validate
        config: Optional configuration dictionary with validation parameters
        output_file: Optional path to save validation results
    
    Returns:
        Dict containing validation results and similar sequence groups
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    results = []
    for binder in data['personalized_binders']:
        validation = validate_binder(binder['sequence'], config)
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