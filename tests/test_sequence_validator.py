"""
Unit tests for sequence validation functionality.
"""

import unittest
import os
import json
from modules.sequence_validator import TherapeuticAntibodyValidator, validate_generated_sequences

class TestSequenceValidator(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.validator = TherapeuticAntibodyValidator()
        
        # Test sequence (based on Adalimumab)
        self.test_sequence = {
            "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS",
            "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIK"
        }
        
    def test_extract_cdrs(self):
        """Test CDR extraction."""
        cdrs = self.validator.extract_cdrs(self.test_sequence['heavy_chain'])
        self.assertIsNotNone(cdrs['cdr1'])
        self.assertIsNotNone(cdrs['cdr2'])
        self.assertIsNotNone(cdrs['cdr3'])
        
    def test_sequence_similarity(self):
        """Test sequence similarity calculation."""
        # Test identical sequences
        similarity = self.validator.calculate_sequence_similarity(
            self.test_sequence['heavy_chain'],
            self.test_sequence['heavy_chain']
        )
        self.assertAlmostEqual(similarity, 1.0, places=2)
        
        # Test different sequences
        different_seq = "X" * len(self.test_sequence['heavy_chain'])
        similarity = self.validator.calculate_sequence_similarity(
            self.test_sequence['heavy_chain'],
            different_seq
        )
        self.assertLess(similarity, 0.5)
        
    def test_find_similar_antibodies(self):
        """Test finding similar antibodies."""
        similar = self.validator.find_similar_antibodies(self.test_sequence)
        self.assertTrue(len(similar) > 0)
        if similar:
            self.assertGreaterEqual(similar[0]['similarity_score'], 0.7)
            
    def test_validate_sequence_metrics(self):
        """Test sequence metrics validation."""
        results = self.validator.validate_sequence_metrics(self.test_sequence)
        self.assertIn('valid', results)
        self.assertIn('warnings', results)
        self.assertIn('metrics', results)
        
    def test_validate_generated_sequences(self):
        """Test batch validation of sequences."""
        sequences = [self.test_sequence]
        results = validate_generated_sequences(sequences)
        self.assertIn('validated_sequences', results)
        self.assertIn('summary', results)
        self.assertEqual(results['summary']['total'], 1)

if __name__ == '__main__':
    unittest.main()