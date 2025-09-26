"""
Unit tests for antibody sequence validation functionality.
"""

import unittest
import os
import json
from modules.antibody_validator import AntibodyValidator, validate_generated_sequences

class TestAntibodyValidator(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.validator = AntibodyValidator()
        
        # Test sequence (based on Adalimumab)
        self.test_sequence = {
            "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS",
            "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIK"
        }
        
    def test_sequence_similarity(self):
        """Test sequence similarity calculation."""
        # Test identical sequences
        similarity = self.validator.calculate_similarity(
            self.test_sequence['heavy_chain'],
            self.test_sequence['heavy_chain']
        )
        self.assertAlmostEqual(similarity, 1.0, places=2)
        
        # Test different sequences
        different_seq = "X" * len(self.test_sequence['heavy_chain'])
        similarity = self.validator.calculate_similarity(
            self.test_sequence['heavy_chain'],
            different_seq
        )
        self.assertLess(similarity, 0.5)
        
    def test_sequence_analysis(self):
        """Test comprehensive sequence analysis."""
        analysis = self.validator.analyze_sequence(self.test_sequence['heavy_chain'])
        
        self.assertIn('length', analysis)
        self.assertIn('molecular_weight', analysis)
        self.assertIn('isoelectric_point', analysis)
        self.assertIn('gravy', analysis)
        self.assertIn('secondary_structure', analysis)
        
        # Test length calculation
        self.assertEqual(analysis['length'], len(self.test_sequence['heavy_chain']))
        
    def test_cdr_validation(self):
        """Test CDR validation."""
        # Test heavy chain CDR3
        cdr3_seq = "AKVSYLSTASSLD"
        validation = self.validator.validate_cdr(cdr3_seq, 'heavy', 3)
        
        self.assertIn('valid', validation)
        self.assertIn('length', validation)
        self.assertIn('properties', validation)
        
    def test_antibody_validation(self):
        """Test full antibody validation."""
        validation = self.validator.validate_antibody(self.test_sequence)
        
        self.assertIn('valid', validation)
        self.assertIn('warnings', validation)
        self.assertIn('analysis', validation)
        
        # Check chain analysis
        self.assertIn('heavy_chain', validation['analysis'])
        self.assertIn('light_chain', validation['analysis'])
        
    def test_find_similar_antibodies(self):
        """Test finding similar antibodies."""
        similar = self.validator.find_similar_antibodies(self.test_sequence)
        
        self.assertIsInstance(similar, list)
        if similar:  # Should find at least one match (itself)
            self.assertIn('name', similar[0])
            self.assertIn('similarity_score', similar[0])
            self.assertGreaterEqual(similar[0]['similarity_score'], 0.7)
            
    def test_batch_validation(self):
        """Test batch validation of sequences."""
        sequences = [self.test_sequence]
        results = validate_generated_sequences(sequences)
        
        self.assertIn('validated_sequences', results)
        self.assertIn('summary', results)
        self.assertEqual(results['summary']['total'], 1)
        
        # Test output file creation
        test_output = 'test_validation_results.json'
        results = validate_generated_sequences(sequences, test_output)
        
        self.assertTrue(os.path.exists(test_output))
        os.remove(test_output)  # Clean up

if __name__ == '__main__':
    unittest.main()