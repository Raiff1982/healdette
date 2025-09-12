"""
Unit tests for sequence validation functions.
"""
import unittest
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from modules.validate_sequences import SequenceValidator

class TestSequenceValidator(unittest.TestCase):
    def setUp(self):
        # Test sequences with known properties
        self.test_sequences = {
            "acidic": "DDDEEEDDEE",  # Should have low pI
            "basic": "KKRRKKKRKR",   # Should have high pI
            "neutral": "GGGAAAGGG",   # Should have middle pI
            "mixed": "DEDKRKDEKR"     # Should have intermediate pI
        }
    
    def test_pi_calculation(self):
        """Test that pI calculations match BioPython's results."""
        for name, seq in self.test_sequences.items():
            validator = SequenceValidator(seq)
            biopython_pi = ProteinAnalysis(seq).isoelectric_point()
            our_pi = validator.calculate_properties()["pI"]
            
            self.assertAlmostEqual(
                our_pi, 
                biopython_pi, 
                places=1,
                msg=f"pI calculation failed for {name} sequence"
            )
    
    def test_gravy_calculation(self):
        """Test GRAVY calculations."""
        hydrophobic = "ILVAAA"  # Should have positive GRAVY
        hydrophilic = "RKDENN"  # Should have negative GRAVY
        
        validator_hydrophobic = SequenceValidator(hydrophobic)
        validator_hydrophilic = SequenceValidator(hydrophilic)
        
        gravy_hydrophobic = validator_hydrophobic.calculate_properties()["GRAVY"]
        gravy_hydrophilic = validator_hydrophilic.calculate_properties()["GRAVY"]
        
        self.assertGreater(gravy_hydrophobic, 0)
        self.assertLess(gravy_hydrophilic, 0)

if __name__ == '__main__':
    unittest.main()