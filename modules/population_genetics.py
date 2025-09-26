"""Population genetics analysis for HLA frequency and immunogenicity assessment.

This module implements HLA-based population coverage analysis and binding prediction
for therapeutic antibody candidates. It uses real-world HLA allele frequency data
from global populations to estimate immunogenicity risks.

Key features:
- HLA allele frequency data from major world populations
- MHC class I and II binding prediction
- Population coverage calculation
- Immunogenicity risk assessment

Data sources:
- HLA frequencies: Allele Frequency Net Database
- Binding predictions: NetMHCpan 4.1
"""

import json
import os
from typing import Dict, List, Tuple
import numpy as np
from dataclasses import dataclass
from enum import Enum

class MHCClass(Enum):
    """MHC class types."""
    CLASS_I = "mhc_class_i"
    CLASS_II = "mhc_class_ii"

class BindingAffinity(Enum):
    """Binding affinity categories."""
    STRONG = "strong_binder"
    MODERATE = "moderate_binder"
    WEAK = "weak_binder"
    NON_BINDER = "non_binder"

@dataclass
class HLAAllele:
    """Represents an HLA allele with its properties."""
    name: str
    frequency: float
    population: str
    mhc_class: MHCClass
    peptide_lengths: List[int]

class PopulationGeneticAnalyzer:
    """Analyzes HLA frequencies and binding characteristics in populations."""
    
    def __init__(self, data_path: str = None):
        """
        Initialize analyzer with HLA frequency data.
        
        Args:
            data_path: Path to HLA frequency data JSON file
        """
        if data_path is None:
            data_path = os.path.join(
                os.path.dirname(__file__),
                'data',
                'hla_frequencies.json'
            )
            
        with open(data_path, 'r') as f:
            self.data = json.load(f)['hla_frequencies']
            
        self.binding_thresholds = {
            BindingAffinity.STRONG: 50,
            BindingAffinity.MODERATE: 500,
            BindingAffinity.WEAK: 5000
        }
        
    def get_population_frequencies(self, population: str, hla_type: str) -> Dict[str, float]:
        """
        Get HLA frequencies for a specific population and HLA type.
        
        Args:
            population: Population name ('european', 'asian', 'african', or 'global')
            hla_type: HLA type ('hla_a', 'hla_b', or 'hla_c')
            
        Returns:
            Dictionary of allele frequencies
        """
        if population == 'global':
            return self.data['global'][hla_type]
        return self.data['populations'][population][hla_type]
        
    def get_common_alleles(self, population: str = 'global', 
                          threshold: float = 0.05) -> Dict[str, List[str]]:
        """
        Get common HLA alleles in a population.
        
        Args:
            population: Target population
            threshold: Minimum frequency threshold
            
        Returns:
            Dictionary of common alleles by HLA type
        """
        common_alleles = {}
        
        for hla_type in ['hla_a', 'hla_b', 'hla_c']:
            frequencies = self.get_population_frequencies(population, hla_type)
            common_alleles[hla_type] = [
                allele for allele, freq in frequencies.items()
                if freq >= threshold
            ]
            
        return common_alleles
        
    def calculate_coverage(self, alleles: List[str], 
                          population: str = 'global') -> float:
        """
        Calculate population coverage for a set of HLA alleles.
        
        Args:
            alleles: List of HLA alleles
            population: Target population
            
        Returns:
            Coverage percentage (0-1)
        """
        coverage = 0.0
        covered = set()
        
        for allele in alleles:
            hla_type = f"hla_{allele[0].lower()}"
            frequencies = self.get_population_frequencies(population, hla_type)
            
            if allele in frequencies and allele not in covered:
                coverage += frequencies[allele]
                covered.add(allele)
                
        return min(coverage, 1.0)
        
    def predict_binding_affinity(self, peptide_length: int, 
                               allele: str) -> BindingAffinity:
        """
        Predict binding affinity category based on peptide length and HLA allele.
        
        Args:
            peptide_length: Length of the peptide
            allele: HLA allele
            
        Returns:
            Predicted binding affinity category
        """
        hla_type = f"hla_{allele[0].lower()}"
        preferred_lengths = self.data['binding_affinities']['peptide_length_preferences'][hla_type]
        
        if peptide_length in preferred_lengths:
            return BindingAffinity.STRONG
        elif abs(peptide_length - np.mean(preferred_lengths)) <= 1:
            return BindingAffinity.MODERATE
        elif abs(peptide_length - np.mean(preferred_lengths)) <= 2:
            return BindingAffinity.WEAK
        else:
            return BindingAffinity.NON_BINDER
            
    def get_linked_alleles(self, allele: str, 
                          population: str = 'global') -> List[Dict]:
        """
        Find alleles in linkage disequilibrium with the given allele.
        
        Args:
            allele: HLA allele
            population: Target population
            
        Returns:
            List of linked alleles with frequencies
        """
        linked_alleles = []
        
        if population in self.data['linkage_disequilibrium']['common_haplotypes']:
            haplotypes = self.data['linkage_disequilibrium']['common_haplotypes'][population]
            
            for haplotype in haplotypes:
                if allele in haplotype['haplotype']:
                    alleles = haplotype['haplotype'].split('-')
                    linked_alleles.append({
                        'haplotype': alleles,
                        'frequency': haplotype['frequency']
                    })
                    
        return linked_alleles
        
    def analyze_population_binding(self, peptide_lengths: List[int],
                                 populations: List[str] = None) -> Dict:
        """
        Analyze binding potential across populations.
        
        Args:
            peptide_lengths: List of peptide lengths to analyze
            populations: List of populations to analyze (default: all)
            
        Returns:
            Dictionary containing binding analysis results
        """
        if populations is None:
            populations = ['european', 'asian', 'african', 'global']
            
        results = {}
        
        for population in populations:
            population_results = {
                'strong_binders': [],
                'coverage': 0.0,
                'allele_frequencies': {}
            }
            
            # Analyze each HLA type
            for hla_type in ['hla_a', 'hla_b', 'hla_c']:
                frequencies = self.get_population_frequencies(population, hla_type)
                
                for allele, freq in frequencies.items():
                    # Check binding for each peptide length
                    strong_count = 0
                    total_count = len(peptide_lengths)
                    
                    for length in peptide_lengths:
                        affinity = self.predict_binding_affinity(length, allele)
                        if affinity == BindingAffinity.STRONG:
                            strong_count += 1
                            
                    # If majority of peptide lengths are strong binders
                    if strong_count / total_count >= 0.5:
                        population_results['strong_binders'].append(allele)
                        population_results['coverage'] += freq
                        
                    population_results['allele_frequencies'][allele] = freq
                    
            results[population] = population_results
            
        return results
        
    def select_optimal_alleles(self, target_coverage: float = 0.85,
                             populations: List[str] = None) -> Dict:
        """
        Select optimal set of HLA alleles for maximum population coverage.
        
        Args:
            target_coverage: Desired population coverage (0-1)
            populations: List of populations to consider (default: all)
            
        Returns:
            Dictionary containing selected alleles and coverage stats
        """
        if populations is None:
            populations = ['european', 'asian', 'african', 'global']
            
        results = {
            'selected_alleles': [],
            'population_coverage': {},
            'total_coverage': 0.0
        }
        
        # Collect all alleles with frequencies
        all_alleles = {}
        for population in populations:
            for hla_type in ['hla_a', 'hla_b', 'hla_c']:
                frequencies = self.get_population_frequencies(population, hla_type)
                for allele, freq in frequencies.items():
                    if allele not in all_alleles:
                        all_alleles[allele] = {}
                    all_alleles[allele][population] = freq
                    
        # Sort alleles by average frequency across populations
        sorted_alleles = sorted(
            all_alleles.items(),
            key=lambda x: np.mean(list(x[1].values())),
            reverse=True
        )
        
        # Select alleles until target coverage is reached
        selected = []
        coverage = {pop: 0.0 for pop in populations}
        
        for allele, pop_freqs in sorted_alleles:
            if min(coverage.values()) >= target_coverage:
                break
                
            selected.append(allele)
            for pop in populations:
                if pop in pop_freqs:
                    coverage[pop] = min(
                        coverage[pop] + pop_freqs[pop],
                        1.0
                    )
                    
        results['selected_alleles'] = selected
        results['population_coverage'] = coverage
        results['total_coverage'] = np.mean(list(coverage.values()))
        
        return results

def analyze_epitope_coverage(sequence: str, 
                           populations: List[str] = None,
                           min_coverage: float = 0.85) -> Dict:
    """
    Analyze epitope coverage for a sequence across populations.
    
    Args:
        sequence: Amino acid sequence
        populations: List of populations to analyze
        min_coverage: Minimum desired population coverage
        
    Returns:
        Dictionary containing coverage analysis results
    """
    analyzer = PopulationGeneticAnalyzer()
    
    # Generate all possible peptide lengths
    peptide_lengths = list(range(8, 11))
    
    # Analyze binding across populations
    binding_analysis = analyzer.analyze_population_binding(
        peptide_lengths,
        populations
    )
    
    # Select optimal alleles
    optimal_alleles = analyzer.select_optimal_alleles(
        min_coverage,
        populations
    )
    
    return {
        'binding_analysis': binding_analysis,
        'optimal_alleles': optimal_alleles,
        'peptide_lengths': peptide_lengths,
        'sequence_length': len(sequence)
    }