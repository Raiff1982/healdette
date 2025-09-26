"""
Unit tests for population genetics functionality.
"""

import unittest
import os
from modules.population_genetics import (
    PopulationGeneticAnalyzer,
    analyze_epitope_coverage,
    MHCClass,
    BindingAffinity
)

class TestPopulationGenetics(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = PopulationGeneticAnalyzer()
        
        # Test sequence (part of Adalimumab CDR3)
        self.test_sequence = "AKVSYLSTASSLD"
        
    def test_population_frequencies(self):
        """Test HLA frequency retrieval."""
        # Test global frequencies
        global_freq = self.analyzer.get_population_frequencies('global', 'hla_a')
        self.assertIsInstance(global_freq, dict)
        self.assertGreater(len(global_freq), 0)
        
        # Test population-specific frequencies
        euro_freq = self.analyzer.get_population_frequencies('european', 'hla_a')
        self.assertIsInstance(euro_freq, dict)
        self.assertGreater(len(euro_freq), 0)
        
    def test_common_alleles(self):
        """Test common allele identification."""
        common = self.analyzer.get_common_alleles(threshold=0.05)
        
        self.assertIn('hla_a', common)
        self.assertIn('hla_b', common)
        self.assertIn('hla_c', common)
        
        # Check that all alleles meet threshold
        for hla_type, alleles in common.items():
            frequencies = self.analyzer.get_population_frequencies('global', hla_type)
            for allele in alleles:
                self.assertGreaterEqual(frequencies[allele], 0.05)
                
    def test_coverage_calculation(self):
        """Test population coverage calculation."""
        alleles = ['A*02:01', 'B*07:02', 'C*07:01']
        coverage = self.analyzer.calculate_coverage(alleles)
        
        self.assertIsInstance(coverage, float)
        self.assertGreaterEqual(coverage, 0.0)
        self.assertLessEqual(coverage, 1.0)
        
    def test_binding_affinity(self):
        """Test binding affinity prediction."""
        # Test optimal length
        affinity = self.analyzer.predict_binding_affinity(9, 'A*02:01')
        self.assertEqual(affinity, BindingAffinity.STRONG)
        
        # Test suboptimal length
        affinity = self.analyzer.predict_binding_affinity(12, 'A*02:01')
        self.assertEqual(affinity, BindingAffinity.NON_BINDER)
        
    def test_linked_alleles(self):
        """Test linkage disequilibrium analysis."""
        linked = self.analyzer.get_linked_alleles('A*01:01', 'european')
        
        self.assertIsInstance(linked, list)
        if linked:  # If there are linked alleles
            self.assertIn('haplotype', linked[0])
            self.assertIn('frequency', linked[0])
            
    def test_population_binding(self):
        """Test population-wide binding analysis."""
        peptide_lengths = [8, 9, 10]
        analysis = self.analyzer.analyze_population_binding(
            peptide_lengths,
            ['european']
        )
        
        self.assertIn('european', analysis)
        self.assertIn('strong_binders', analysis['european'])
        self.assertIn('coverage', analysis['european'])
        
    def test_optimal_allele_selection(self):
        """Test optimal allele selection."""
        results = self.analyzer.select_optimal_alleles(
            target_coverage=0.85,
            populations=['european', 'asian']
        )
        
        self.assertIn('selected_alleles', results)
        self.assertIn('population_coverage', results)
        self.assertIn('total_coverage', results)
        
        # Check coverage values
        self.assertGreaterEqual(results['total_coverage'], 0.0)
        self.assertLessEqual(results['total_coverage'], 1.0)
        
    def test_epitope_coverage(self):
        """Test epitope coverage analysis."""
        results = analyze_epitope_coverage(
            self.test_sequence,
            populations=['european', 'asian'],
            min_coverage=0.85
        )
        
        self.assertIn('binding_analysis', results)
        self.assertIn('optimal_alleles', results)
        self.assertIn('peptide_lengths', results)
        self.assertEqual(results['sequence_length'], len(self.test_sequence))

if __name__ == '__main__':
    unittest.main()