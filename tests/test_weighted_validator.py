"""
Unit tests for the weighted sequence validator.
"""

import unittest
import json
import os
from typing import Dict
from modules.weighted_validator import WeightedSequenceValidator

class TestWeightedValidator(unittest.TestCase):
    def setUp(self):
        """Set up test cases."""
        self.test_sequence = "WYRKFGKHWFRYKHFWYRFKHWFYKHFRWYKFHWYKHFWYRFKHWYRKHFW"
        self.test_config = {
            "global_params": {
                "sequence_length": {
                    "min": 40,
                    "max": 70
                }
            },
            "populations": {
                "french_german": {
                    "ancestry_weight": 0.6,
                    "binding_motifs": ["WY", "RF", "KH", "YF"],
                    "biophysical_params": {
                        "aromatic_content": {
                            "min": 16,
                            "max": 28
                        },
                        "hydrophobic_content": {
                            "min": 33,
                            "max": 43
                        }
                    }
                },
                "finnish": {
                    "ancestry_weight": 0.4,
                    "binding_motifs": ["WH", "RF", "KY", "FF"],
                    "biophysical_params": {
                        "aromatic_content": {
                            "min": 14,
                            "max": 26
                        },
                        "hydrophobic_content": {
                            "min": 32,
                            "max": 42
                        }
                    }
                }
            }
        }
        
    def test_parameter_blending(self):
        """Test that parameters are correctly blended based on ancestry weights."""
        validator = WeightedSequenceValidator(self.test_sequence, self.test_config)
        params = validator._blend_parameters()
        
        # Expected values based on weights (0.6 and 0.4)
        expected_aromatic_min = round(16 * 0.6 + 14 * 0.4, 2)
        expected_aromatic_max = round(28 * 0.6 + 26 * 0.4, 2)
        
        self.assertEqual(params['aromatic_content']['min'], expected_aromatic_min)
        self.assertEqual(params['aromatic_content']['max'], expected_aromatic_max)
        
    def test_binding_motifs(self):
        """Test binding motif detection and scoring."""
        validator = WeightedSequenceValidator(self.test_sequence, self.test_config)
        motif_results = validator._check_binding_motifs()
        
        self.assertTrue('scores' in motif_results)
        self.assertTrue('total_score' in motif_results)
        self.assertTrue(0 <= motif_results['total_score'] <= 1)
        
    def test_population_scoring(self):
        """Test population-specific scoring."""
        validator = WeightedSequenceValidator(self.test_sequence, self.test_config)
        results = validator.validate_sequence()
        
        self.assertTrue('population_scores' in results)
        self.assertTrue('french_german' in results['population_scores'])
        self.assertTrue('finnish' in results['population_scores'])
        
        for pop_scores in results['population_scores'].values():
            self.assertTrue(0 <= pop_scores['score'] <= 1)
            
    def test_sequence_validation(self):
        """Test complete sequence validation."""
        validator = WeightedSequenceValidator(self.test_sequence, self.test_config)
        results = validator.validate_sequence()
        
        self.assertTrue('valid' in results)
        self.assertTrue('warnings' in results)
        self.assertTrue('metrics' in results)
        self.assertTrue('population_scores' in results)
        
        # Check metric calculations
        self.assertTrue('aromatic_content' in results['metrics'])
        self.assertTrue('hydrophobic_content' in results['metrics'])
        self.assertTrue('binding_motifs' in results['metrics'])

if __name__ == '__main__':
    unittest.main()