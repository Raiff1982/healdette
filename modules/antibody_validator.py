"""
Enhanced sequence validation module for antibody therapeutics.
"""

import json
import os
from typing import Dict, List, Tuple
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class AntibodyValidator:
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
        self.config = {
            "sequence_properties": {
                "min_length": 100,
                "max_length": 500,
                "max_hydrophobicity": 0.5,
                "max_instability": 40.0
            },
            "cdrs": {
                "heavy": {
                    "cdr1_length": [5, 10],
                    "cdr2_length": [12, 17],
                    "cdr3_length": [9, 17]
                },
                "light": {
                    "cdr1_length": [10, 12],
                    "cdr2_length": [7, 8],
                    "cdr3_length": [8, 10]
                }
            }
        }

    def calculate_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence similarity using BioPython's implementation.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
        
        Returns:
            Similarity score between 0 and 1
        """
        # Define amino acid similarity groups
        similar_groups = [
            set('ILVM'),   # Aliphatic
            set('FYW'),    # Aromatic
            set('KRH'),    # Basic
            set('DE'),     # Acidic
            set('STNQ'),   # Polar
            set('AG'),     # Small
            set('C'),      # Cysteine
            set('P')       # Proline
        ]
        
        def aa_similarity(aa1: str, aa2: str) -> float:
            if aa1 == aa2:
                return 1.0
            for group in similar_groups:
                if aa1 in group and aa2 in group:
                    return 0.5
            return 0.0
        
        min_length = min(len(seq1), len(seq2))
        scores = [aa_similarity(a, b) for a, b in zip(seq1[:min_length], seq2[:min_length])]
        length_penalty = abs(len(seq1) - len(seq2)) * 0.1
        
        return max(0.0, min(1.0, sum(scores) / len(scores) - length_penalty)) if scores else 0.0

    def analyze_sequence(self, sequence: str) -> Dict:
        """
        Perform comprehensive sequence analysis using BioPython.
        
        Args:
            sequence: Amino acid sequence to analyze
            
        Returns:
            Dictionary containing analysis results
        """
        try:
            analyzer = ProteinAnalysis(sequence)
            
            return {
                "length": len(sequence),
                "molecular_weight": round(analyzer.molecular_weight(), 2),
                "aromaticity": round(analyzer.aromaticity(), 3),
                "instability_index": round(analyzer.instability_index(), 2),
                "isoelectric_point": round(analyzer.isoelectric_point(), 2),
                "gravy": round(analyzer.gravy(), 3),
                "secondary_structure": {
                    "helix": round(analyzer.secondary_structure_fraction()[0], 3),
                    "sheet": round(analyzer.secondary_structure_fraction()[1], 3),
                    "coil": round(analyzer.secondary_structure_fraction()[2], 3)
                },
                "amino_acid_composition": {
                    aa: round(percent, 3)
                    for aa, percent in analyzer.get_amino_acids_percent().items()
                }
            }
        except Exception as e:
            return {
                "error": f"Analysis failed: {str(e)}",
                "length": len(sequence)
            }

    def validate_cdr(self, cdr_seq: str, chain_type: str, cdr_num: int) -> Dict:
        """
        Validate a CDR sequence against known patterns.
        
        Args:
            cdr_seq: CDR sequence to validate
            chain_type: 'heavy' or 'light'
            cdr_num: CDR number (1, 2, or 3)
            
        Returns:
            Validation results
        """
        length_range = self.config['cdrs'][chain_type][f'cdr{cdr_num}_length']
        
        results = {
            "valid": True,
            "length": len(cdr_seq),
            "warnings": []
        }
        
        # Length validation
        if len(cdr_seq) < length_range[0] or len(cdr_seq) > length_range[1]:
            results["valid"] = False
            results["warnings"].append(
                f"Length {len(cdr_seq)} outside valid range {length_range}"
            )
        
        # Analyze composition
        try:
            analyzer = ProteinAnalysis(cdr_seq)
            results["properties"] = {
                "hydrophobicity": round(analyzer.gravy(), 3),
                "isoelectric_point": round(analyzer.isoelectric_point(), 2),
                "aromatic_fraction": round(sum(
                    cdr_seq.count(aa) for aa in 'FWY'
                ) / len(cdr_seq), 3)
            }
            
            # CDR-specific checks
            if cdr_num == 3:  # CDR3-specific validation
                if results["properties"]["aromatic_fraction"] < 0.1:
                    results["warnings"].append("Low aromatic content in CDR3")
            
        except Exception:
            results["properties"] = None
            results["warnings"].append("Could not analyze CDR properties")
        
        return results

    def validate_antibody(self, sequence: Dict[str, str]) -> Dict:
        """
        Perform comprehensive antibody validation.
        
        Args:
            sequence: Dictionary with 'heavy_chain' and 'light_chain' sequences
            
        Returns:
            Validation results
        """
        results = {
            "valid": True,
            "warnings": [],
            "analysis": {}
        }
        
        # Validate each chain
        for chain in ['heavy_chain', 'light_chain']:
            if chain not in sequence:
                results["valid"] = False
                results["warnings"].append(f"Missing {chain}")
                continue
                
            # Analyze sequence
            chain_analysis = self.analyze_sequence(sequence[chain])
            results["analysis"][chain] = chain_analysis
            
            # Length validation
            length_range = self.metrics["sequence_length"][chain]
            if chain_analysis["length"] < length_range["min"] or chain_analysis["length"] > length_range["max"]:
                results["valid"] = False
                results["warnings"].append(
                    f"{chain} length {chain_analysis['length']} outside valid range "
                    f"({length_range['min']}-{length_range['max']})"
                )
        
        # Find similar antibodies
        if results["valid"]:
            similar = self.find_similar_antibodies(sequence)
            if similar:
                results["similar_antibodies"] = similar
                
                # Warning if too similar to existing antibody
                if similar[0]["similarity_score"] > 0.9:
                    results["warnings"].append(
                        f"Very high similarity ({similar[0]['similarity_score']:.3f}) "
                        f"to existing antibody {similar[0]['name']}"
                    )
        
        return results

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
            heavy_similarity = self.calculate_similarity(
                sequence['heavy_chain'],
                antibody['sequence']['heavy_chain']
            )
            
            light_similarity = self.calculate_similarity(
                sequence['light_chain'],
                antibody['sequence']['light_chain']
            )
            
            # Weight heavy chain more in similarity calculation
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
    validator = AntibodyValidator()
    
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
        validation = validator.validate_antibody(sequence)
        
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