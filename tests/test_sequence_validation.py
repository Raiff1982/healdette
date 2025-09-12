"""
Unit tests for sequence validation module, including pI calculation and complexity analysis.
"""

import unittest
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sys
import os

# Add the parent directory to the path so we can import our modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.validate_sequences import SequenceValidator, validate_binder

class TestSequenceValidation(unittest.TestCase):
    def setUp(self):
        # Test sequences for different validation aspects
        self.test_sequences = {
            'acidic': 'DDDEEEDDDEEE',  # Should have low pI
            'basic': 'KKRRKKKRRKRK',   # Should have high pI
            'neutral': 'GGGGGGGGGGGG', # Should be near neutral
            'mixed': 'KDKEFGYWAPTS',   # Mix of amino acids
            'real_peptide': 'MKKSFWLVLLVALNLWIKANA',  # Realistic signal peptide
            # Sequences for complexity testing
            'homopolymer': 'MKAAAAATWLVLLVALNLWIKANA',  # Has AAAAA run
            'aqp_heavy': 'MAPQAPQAPQAPQAPQAPQTWLVL',    # High A/Q/P content
            'low_complexity': 'AQPAQPAQPAQPAQPAQPAQP',   # Very repetitive
            'good_binder': 'MKKSFWLVLLCALNLWIKANACR',    # Well-balanced sequence
            # Sequences for cysteine pattern testing
            'terminal_pair': 'MKKCFWLVLLVALNLWIKANACT',   # Terminal cysteine pair
            'ladder_motif': 'MKCDEFCGHICKLMCNOPQCR',     # Evenly spaced cysteines
            'odd_cysteines': 'MKCDEFCGHICKLMNOPQRS',     # Odd number of cysteines
            'optimal_scaffold': 'MKCDEFGHICKLMNOPQCRSTC' # Good scaffold pattern (4 cysteines)
        }
        
    def test_pi_calculation_range(self):
        """
        Test pI calculation for Codette binder requirements:
        - Acidic sequences (pI < 5): Important for stability in physiological conditions
        - Neutral sequences (6 < pI < 8): Typical for well-behaved binders
        - Basic sequences (pI > 9): Important for target binding
        """
        test_cases = [
            ('DDDEEEDDDEEE', 'acidic', lambda x: x < 5),
            ('GGGGGGGGGGGG', 'neutral', lambda x: 5 <= x <= 8),  # Broader neutral range for Codette binders
            ('KKRRKKKRRKRK', 'basic', lambda x: x > 9)
        ]
        
        for seq, category, validator_func in test_cases:
            validator = SequenceValidator(seq)
            pi = validator.calculate_properties()['pI']
            
            self.assertTrue(
                validator_func(pi),
                f"pI {pi} for {category} sequence {seq} is outside expected range"
            )
            
    def test_charge_ph_relationships(self):
        """Test specific charge/pH relationships required for Codette binders"""
        # Test acidic sequence
        acidic_seq = 'DDDEEEDDDEEE'
        validator = SequenceValidator(acidic_seq)
        
        # At pH 7.4 (physiological), acidic sequences should have significant negative charge
        charge_phys = validator.charge_at_ph(7.4)
        self.assertLess(
            charge_phys, 
            -5.0,
            f"Acidic sequence charge at pH 7.4 ({charge_phys}) not negative enough"
        )
        
        # Basic sequence at physiological pH
        basic_seq = 'KKRRKKKRRKRK'
        validator = SequenceValidator(basic_seq)
        charge_phys = validator.charge_at_ph(7.4)
        self.assertGreater(
            charge_phys,
            5.0,
            f"Basic sequence charge at pH 7.4 ({charge_phys}) not positive enough"
        )
    
    def test_sequence_complexity(self):
        """Test sequence complexity analysis"""
        # Test homopolymer detection
        validator = SequenceValidator(self.test_sequences['homopolymer'])
        complexity = validator.analyze_complexity()
        self.assertTrue(
            any(run['amino_acid'] == 'A' and run['length'] >= 5 
                for run in complexity['homopolymer_runs']),
            "Failed to detect AAAAA homopolymer run"
        )
        
        # Test A/Q/P-heavy regions
        validator = SequenceValidator(self.test_sequences['aqp_heavy'])
        complexity = validator.analyze_complexity()
        self.assertTrue(
            complexity['warnings']['high_aqp'],
            "Failed to detect high A/Q/P content"
        )
        self.assertGreater(
            len(complexity['aqp_heavy_regions']),
            0,
            "Failed to identify A/Q/P-heavy regions"
        )
        
        # Test sequence entropy
        validator = SequenceValidator(self.test_sequences['low_complexity'])
        complexity = validator.analyze_complexity()
        self.assertTrue(
            complexity['warnings']['low_complexity'],
            "Failed to detect low complexity sequence"
        )
        self.assertLess(
            complexity['sequence_entropy'],
            3.0,
            "Low complexity sequence has unexpectedly high entropy"
        )
        
        # Test well-balanced sequence
        validator = SequenceValidator(self.test_sequences['good_binder'])
        complexity = validator.analyze_complexity()
        self.assertFalse(
            any(complexity['warnings'].values()),
            "Good binder sequence incorrectly flagged with warnings"
        )
        
    def test_cysteine_analysis(self):
        """Test enhanced cysteine pattern analysis for binder scaffolds"""
        # Test terminal pair pattern
        validator = SequenceValidator(self.test_sequences['terminal_pair'])
        analysis = validator.analyze_cysteines()
        self.assertTrue(
            analysis['patterns']['motifs']['terminal_pair'],
            "Failed to detect terminal cysteine pair pattern"
        )
        self.assertTrue(
            analysis['scaffold_evaluation']['suitable_scaffold'],
            "Terminal pair pattern not recognized as suitable scaffold"
        )
        
        # Test ladder motif
        validator = SequenceValidator(self.test_sequences['ladder_motif'])
        analysis = validator.analyze_cysteines()
        self.assertTrue(
            analysis['patterns']['motifs']['ladder'],
            "Failed to detect ladder-like cysteine pattern"
        )
        self.assertTrue(
            analysis['scaffold_evaluation']['suitable_scaffold'],
            "Ladder pattern not recognized as suitable scaffold"
        )
        
        # Test odd number of cysteines
        validator = SequenceValidator(self.test_sequences['odd_cysteines'])
        analysis = validator.analyze_cysteines()
        self.assertFalse(
            analysis['patterns']['paired'],
            "Odd number of cysteines incorrectly marked as paired"
        )
        self.assertTrue(
            any("Odd number of cysteines" in warning 
                for warning in analysis['warnings'] if warning),
            "No warning for odd number of cysteines"
        )
        
        # Test optimal scaffold
        test_seq = self.test_sequences['optimal_scaffold']
        cys_count = test_seq.count('C')
        print(f"\nDebug - Optimal scaffold sequence: {test_seq}")
        print(f"Debug - Cysteine count: {cys_count}")
        print(f"Debug - Cysteine positions: {[i for i, aa in enumerate(test_seq) if aa == 'C']}")
        
        validator = SequenceValidator(self.test_sequences['optimal_scaffold'])
        analysis = validator.analyze_cysteines()
        
        print(f"Debug - Analysis result: {analysis}")
        
        self.assertTrue(
            analysis['scaffold_evaluation']['optimal_count'],
            "Optimal cysteine count not recognized"
        )
        self.assertTrue(
            analysis['scaffold_evaluation']['well_distributed'],
            "Well-distributed cysteines not recognized"
        )
        self.assertTrue(
            analysis['scaffold_evaluation']['suitable_scaffold'],
            "Optimal scaffold pattern not recognized"
        )
        warnings = [w for w in analysis['warnings'] if w]
        self.assertEqual(
            len(warnings),
            0,
            f"Unexpected warnings for optimal scaffold: {warnings}"
        )

    def test_comprehensive_validation(self):
        """Test the complete binder validation process"""
        # Test problematic sequence
        bad_sequence = 'AAAAAQQQQPPPPPAAAQQP'
        result = validate_binder(bad_sequence)
        
        self.assertTrue('warnings' in result)
        self.assertGreater(len(result['warnings']), 0)
        self.assertFalse(result['is_valid'])
        
        # Warnings should mention specific issues
        warning_text = ' '.join(result['warnings']).lower()
        self.assertTrue(
            any(term in warning_text for term in ['homopolymer', 'low complexity', 'high a/q/p']),
            "Warnings don't specify sequence complexity issues"
        )
        
        # Test good sequence
        good_sequence = self.test_sequences['good_binder']
        result = validate_binder(good_sequence)
        
        self.assertTrue('warnings' in result)
        self.assertEqual(len(result['warnings']), 0)
        self.assertTrue(result['is_valid'])
            
    def test_charge_calculation(self):
        """Test charge calculation at specific pH values."""
        test_cases = [
            ('DDDEEEDDDEEE', 7.0, -12),  # Acidic sequence at neutral pH
            ('KKRRKKKRRKRK', 7.0, 12),   # Basic sequence at neutral pH
            ('GGGGGGGGGGGG', 7.0, 0),    # Neutral sequence
        ]
        
        for seq, ph, expected_charge in test_cases:
            validator = SequenceValidator(seq)
            charge = validator.charge_at_ph(ph)
            self.assertAlmostEqual(
                charge, 
                expected_charge, 
                places=0,
                msg=f"Charge calculation incorrect for sequence {seq} at pH {ph}"
            )

if __name__ == '__main__':
    unittest.main()