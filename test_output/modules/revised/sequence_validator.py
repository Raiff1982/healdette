"""Sequence validation module for antibody generation."""

from typing import Dict, List
from collections import Counter
import re

class SequenceValidator:
    """Validates antibody sequences for quality and structural properties."""
    
    def __init__(self, config: Dict = None):
        """Initialize validator with configurable parameters."""
        default_config = {
            'max_aa_frequency': 0.35,  # Maximum frequency of any single amino acid
            'min_unique_aa': 5,  # Minimum number of unique amino acids
            'max_homopolymer': 4,  # Maximum consecutive repeat of same amino acid
            'min_hydrophobic': 0.15,  # Minimum hydrophobic content
            'max_hydrophobic': 0.65,  # Maximum hydrophobic content
            'charged_ratio_range': (0.1, 0.4),  # Min/max ratio of charged residues
        }
        self.config = config or default_config
        self.hydrophobic = set('AILMFWYV')
        self.charged = set('DEKRH')
        
    def analyze_sequence(self, sequence: str) -> Dict:
        """Analyze sequence properties and validate against criteria."""
        if not sequence or not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence):
            return {
                'passes_validation': False,
                'error': 'Invalid sequence - contains non-standard amino acids'
            }
        
        analysis = {
            'length': len(sequence),
            'aa_composition': self._get_aa_composition(sequence),
            'homopolymer_runs': self._find_homopolymer_runs(sequence),
            'hydrophobic_content': self._calculate_hydrophobic_content(sequence),
            'charge_properties': self._analyze_charge_properties(sequence),
            'structure_metrics': self._analyze_structure(sequence)
        }
        
        # Validate against criteria
        validation = {
            'max_aa_freq_ok': max(analysis['aa_composition'].values()) <= self.config['max_aa_frequency'],
            'unique_aa_ok': len(analysis['aa_composition']) >= self.config['min_unique_aa'],
            'homopolymer_ok': all(len(run) <= self.config['max_homopolymer'] for run in analysis['homopolymer_runs']),
            'hydrophobic_ok': (self.config['min_hydrophobic'] <= analysis['hydrophobic_content'] <= 
                             self.config['max_hydrophobic']),
            'charge_ok': (self.config['charged_ratio_range'][0] <= 
                         analysis['charge_properties']['charged_ratio'] <= 
                         self.config['charged_ratio_range'][1])
        }
        
        analysis['validation'] = validation
        analysis['passes_validation'] = all(validation.values())
        
        return analysis
    
    def _get_aa_composition(self, sequence: str) -> Dict[str, float]:
        """Calculate amino acid composition as frequencies."""
        counts = Counter(sequence)
        length = len(sequence)
        return {aa: count/length for aa, count in counts.items()}
    
    def _find_homopolymer_runs(self, sequence: str) -> List[str]:
        """Find runs of identical amino acids."""
        pattern = r'(.)\1{' + str(self.config['max_homopolymer']-1) + ',}'
        return [match.group() for match in re.finditer(pattern, sequence)]
    
    def _calculate_hydrophobic_content(self, sequence: str) -> float:
        """Calculate fraction of hydrophobic residues."""
        return sum(aa in self.hydrophobic for aa in sequence) / len(sequence)
    
    def _analyze_charge_properties(self, sequence: str) -> Dict:
        """Analyze charge distribution and properties."""
        charged_count = sum(aa in self.charged for aa in sequence)
        return {
            'charged_ratio': charged_count / len(sequence),
            'net_charge': sum((aa in 'KR') - (aa in 'DE') for aa in sequence)
        }
    
    def _analyze_structure(self, sequence: str) -> Dict:
        """Calculate basic structural properties."""
        # Simplified secondary structure propensity
        helix_prone = set('AMLKRH')
        sheet_prone = set('VIFYW')
        
        return {
            'helix_propensity': sum(aa in helix_prone for aa in sequence) / len(sequence),
            'sheet_propensity': sum(aa in sheet_prone for aa in sequence) / len(sequence)
        }