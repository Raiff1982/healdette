"""
Enhanced sequence validator with ancestry weight-based parameter blending.
"""

import re
import json
import math
import datetime
from typing import Dict, List, Tuple, Optional
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class WeightedSequenceValidator:
    def __init__(self, sequence: str, population_config: Dict):
        """
        Initialize validator with population-specific parameters.
        
        Args:
            sequence: The amino acid sequence to validate
            population_config: Configuration with population-specific parameters and weights
        """
        self.sequence = sequence.upper()
        self.populations = population_config.get('populations', {})
        self.global_params = population_config.get('global_params', {})
        
    def _blend_parameters(self) -> Dict:
        """
        Blend parameters from different populations based on ancestry weights.
        """
        blended_params = {
            'aromatic_content': {'min': 0, 'max': 0},
            'hydrophobic_content': {'min': 0, 'max': 0},
            'net_charge': {'min': 0, 'max': 0}
        }
        
        total_weight = sum(pop.get('ancestry_weight', 0) for pop in self.populations.values())
        
        if total_weight == 0:
            return None
        
        # Blend parameters based on ancestry weights
        for pop_name, pop_data in self.populations.items():
            weight = pop_data.get('ancestry_weight', 0) / total_weight
            params = pop_data.get('biophysical_params', {})
            
            for param in blended_params:
                if param in params:
                    blended_params[param]['min'] += params[param]['min'] * weight
                    blended_params[param]['max'] += params[param]['max'] * weight
        
        # Round the blended parameters
        for param in blended_params:
            blended_params[param]['min'] = round(blended_params[param]['min'], 2)
            blended_params[param]['max'] = round(blended_params[param]['max'], 2)
            
        return blended_params
    
    def _check_binding_motifs(self) -> Dict:
        """
        Check for population-specific binding motifs with weighted importance.
        """
        motif_scores = {}
        total_weight = sum(pop.get('ancestry_weight', 0) for pop in self.populations.values())
        
        if total_weight == 0:
            return {'valid': False, 'message': 'No population weights defined'}
            
        for pop_name, pop_data in self.populations.items():
            weight = pop_data.get('ancestry_weight', 0) / total_weight
            motifs = pop_data.get('binding_motifs', [])
            
            pop_score = 0
            for motif in motifs:
                if motif in self.sequence:
                    pop_score += 1
            
            motif_scores[pop_name] = {
                'score': pop_score / len(motifs) if motifs else 0,
                'weighted_score': (pop_score / len(motifs) if motifs else 0) * weight
            }
        
        total_score = sum(score['weighted_score'] for score in motif_scores.values())
        
        return {
            'valid': total_score >= 0.3,  # At least 30% weighted motif match
            'scores': motif_scores,
            'total_score': total_score
        }
    
    def validate_sequence(self) -> Dict:
        """
        Perform ancestry-weighted validation of the sequence.
        """
        results = {
            'valid': True,
            'warnings': [],
            'metrics': {},
            'population_scores': {}
        }
        
        # Get blended parameters
        params = self._blend_parameters()
        if not params:
            results['valid'] = False
            results['warnings'].append('No valid population parameters found')
            return results
            
        # Basic sequence validation
        if not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in self.sequence):
            results['valid'] = False
            results['warnings'].append('Invalid amino acids present')
            return results
            
        # Length validation
        length = len(self.sequence)
        if length < self.global_params.get('sequence_length', {}).get('min', 40) or \
           length > self.global_params.get('sequence_length', {}).get('max', 70):
            results['valid'] = False
            results['warnings'].append(f'Length {length} outside allowed range')
        
        # Analyze with BioPython's ProteinAnalysis
        analysis = ProteinAnalysis(self.sequence)
        
        # Check aromatic content
        aromatic_aas = "FWY"
        aromatic_content = sum(self.sequence.count(aa) for aa in aromatic_aas) / length
        results['metrics']['aromatic_content'] = round(aromatic_content * 100, 2)
        
        if not (params['aromatic_content']['min']/100 <= aromatic_content <= params['aromatic_content']['max']/100):
            results['warnings'].append('Aromatic content outside blended range')
            results['valid'] = False
            
        # Check hydrophobic content
        hydrophobic_aas = "AILMFWYV"
        hydrophobic_content = sum(self.sequence.count(aa) for aa in hydrophobic_aas) / length
        results['metrics']['hydrophobic_content'] = round(hydrophobic_content * 100, 2)
        
        if not (params['hydrophobic_content']['min']/100 <= hydrophobic_content <= params['hydrophobic_content']['max']/100):
            results['warnings'].append('Hydrophobic content outside blended range')
            results['valid'] = False
            
        # Check binding motifs
        motif_results = self._check_binding_motifs()
        results['metrics']['binding_motifs'] = motif_results
        
        if not motif_results['valid']:
            results['warnings'].append('Insufficient population-specific binding motifs')
            results['valid'] = False
            
        # Add population-specific scores
        for pop_name, pop_data in self.populations.items():
            pop_score = self._calculate_population_score(pop_data)
            results['population_scores'][pop_name] = {
                'score': pop_score,
                'weight': pop_data.get('ancestry_weight', 0)
            }
            
        return results
    
    def _calculate_population_score(self, pop_data: Dict) -> float:
        """
        Calculate how well the sequence matches a specific population's requirements.
        """
        score = 0
        params = pop_data.get('biophysical_params', {})
        
        # Check aromatic content
        aromatic_content = sum(self.sequence.count(aa) for aa in "FWY") / len(self.sequence) * 100
        if params.get('aromatic_content'):
            if params['aromatic_content']['min'] <= aromatic_content <= params['aromatic_content']['max']:
                score += 0.3
                
        # Check hydrophobic content
        hydrophobic_content = sum(self.sequence.count(aa) for aa in "AILMFWYV") / len(self.sequence) * 100
        if params.get('hydrophobic_content'):
            if params['hydrophobic_content']['min'] <= hydrophobic_content <= params['hydrophobic_content']['max']:
                score += 0.3
                
        # Check binding motifs
        motifs = pop_data.get('binding_motifs', [])
        if motifs:
            motif_matches = sum(1 for motif in motifs if motif in self.sequence)
            score += 0.4 * (motif_matches / len(motifs))
            
        return round(score, 3)

def validate_sequence_with_ancestry(sequence: str, config_path: str) -> Dict:
    """
    Validate a sequence using ancestry-weighted parameters.
    
    Args:
        sequence: The amino acid sequence to validate
        config_path: Path to the configuration file with population parameters
        
    Returns:
        Dict containing validation results
    """
    with open(config_path, 'r') as f:
        config = json.load(f)
        
    validator = WeightedSequenceValidator(sequence, config)
    results = validator.validate_sequence()
    
    # Add metadata
    results['timestamp'] = datetime.datetime.now().isoformat()
    results['validator_version'] = '3.0.0'
    
    return results