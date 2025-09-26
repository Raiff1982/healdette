"""Sequence validation module for antibody generation."""

from typing import Dict, List
from collections import Counter
import re

class SequenceValidator:
    """Validates antibody sequences for quality and structural properties."""
    
    def __init__(self, config: Dict = None):
        """Initialize validator with configurable parameters."""
        default_config = {
            'max_aa_frequency': 0.4,  # Maximum frequency of any single amino acid
            'min_unique_aa': 4,  # Minimum number of unique amino acids
            'max_homopolymer': 5,  # Maximum consecutive repeat of same amino acid
            'min_hydrophobic': 0.1,  # Minimum hydrophobic content
            'max_hydrophobic': 0.7,  # Maximum hydrophobic content
            'charged_ratio_range': (0.05, 0.5),  # Min/max ratio of charged residues
            'min_length': 8,  # Minimum length for HLA binding
            'max_length': 100,  # Maximum length for stability
            'min_aromatic': 0.05,  # Minimum aromatic content for HLA binding
            'max_aromatic': 0.3,  # Maximum aromatic content
            'celtic_motifs': ['YF', 'WY', 'RF', 'KW']  # Common Celtic HLA-binding motifs
        }
        self.config = config or default_config
        self.hydrophobic = set('AILMFWYV')
        self.charged = set('DEKRH')
        
    def analyze_sequence(self, sequence: str) -> Dict:
        """Analyze sequence properties and validate against criteria."""
        print(f"\nAnalyzing sequence: {sequence}")
        
        if not sequence or not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence):
            print("  Invalid sequence - contains non-standard amino acids")
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
            'structure_metrics': self._analyze_structure(sequence),
            'aromatic_content': self._calculate_aromatic_content(sequence),
            'celtic_binding_motifs': self._find_celtic_motifs(sequence)
        }
        
        print("  Analysis results:")
        print(f"    Length: {analysis['length']}")
        print(f"    AA composition: {analysis['aa_composition']}")
        print(f"    Homopolymer runs: {analysis['homopolymer_runs']}")
        print(f"    Hydrophobic content: {analysis['hydrophobic_content']:.2f}")
        print(f"    Aromatic content: {analysis['aromatic_content']:.2f}")
        print(f"    Celtic binding motifs: {analysis['celtic_binding_motifs']}")
        print(f"    Charge properties: {analysis['charge_properties']}")
        print(f"    Structure metrics: {analysis['structure_metrics']}")
        
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
        
    def _calculate_aromatic_content(self, sequence: str) -> float:
        """Calculate fraction of aromatic residues."""
        aromatics = set('FWY')
        return sum(aa in aromatics for aa in sequence) / len(sequence)
        
    def _find_celtic_motifs(self, sequence: str) -> List[str]:
        """Find Celtic HLA-binding motifs in sequence."""
        found_motifs = []
        for motif in self.config['celtic_motifs']:
            if motif in sequence:
                found_motifs.append(motif)
        return found_motifs