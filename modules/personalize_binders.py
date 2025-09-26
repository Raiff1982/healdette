
import json
import os
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PersonalizationEngine:
    def __init__(self):
        self.hla_data = self._load_hla_data()
        self.celtic_motifs = {
            'RF': 0.15,  # Celtic-specific binding motif weights
            'KW': 0.12,
            'WY': 0.10,
            'YF': 0.08
        }
        
    def _load_hla_data(self):
        """Load HLA frequency data from JSON file."""
        try:
            data_path = Path(__file__).parent.parent / 'data' / 'hla_frequencies.json'
            with open(data_path) as f:
                return json.load(f)['hla_frequencies']
        except Exception as e:
            logger.warning(f"Failed to load HLA data: {e}")
            return None

    def calculate_celtic_score(self, sequence):
        """Calculate Celtic-specific binding score."""
        score = 0
        for motif, weight in self.celtic_motifs.items():
            if motif in sequence:
                score += weight
        return score

    def calculate_population_coverage(self, sequence, ancestry_profile):
        """Calculate population-specific coverage."""
        if not self.hla_data or not self.hla_data['populations']:
            return 0.0
            
        coverage = 0.0
        population_weights = {
            'celtic_british': 0.6,
            'irish_dominant': 0.25,
            'scottish': 0.15
        }
        
        for pop, weight in population_weights.items():
            if pop in self.hla_data['populations']:
                pop_data = self.hla_data['populations'][pop]
                # Calculate weighted coverage for each HLA type
                for hla_type in ['hla_a', 'hla_b', 'hla_c']:
                    if hla_type in pop_data:
                        coverage += sum(pop_data[hla_type].values()) * weight
                        
        return min(coverage, 1.0)

    def personalize_sequence(self, sequence, patient_data):
        """
        Personalize antibody sequence for Celtic/British ancestry.
        
        Args:
            sequence: Original antibody sequence
            patient_data: Dict containing patient's immunogenetic profile
        
        Returns:
            Dict with personalization results
        """
        celtic_score = self.calculate_celtic_score(sequence)
        population_coverage = self.calculate_population_coverage(
            sequence, 
            patient_data.get('ancestry_profile', ['celtic_british'])
        )
        
        # Calculate metabolic adjustment
        metabolic_rate = float(patient_data.get('metabolic_rate', 1.0))
        metabolic_factor = 1.0 / metabolic_rate
        
        # Calculate final score
        base_score = 0.4
        final_score = (
            base_score +
            (celtic_score * 0.3) +
            (population_coverage * 0.2) +
            (metabolic_factor * 0.1)
        )
        
        return {
            'sequence': sequence,
            'celtic_score': round(celtic_score, 3),
            'population_coverage': round(population_coverage, 3),
            'metabolic_factor': round(metabolic_factor, 3),
            'final_score': round(final_score, 3),
            'validation': {
                'celtic_motifs_found': [m for m in self.celtic_motifs if m in sequence],
                'coverage_level': 'high' if population_coverage > 0.7 else 'moderate',
                'metabolic_adjustment': 'standard' if 0.9 <= metabolic_rate <= 1.1 else 'custom'
            }
        }

def personalize_binders(validated_input, patient_data):
    """
    Main function to personalize antibody candidates.
    
    Args:
        validated_input: Dict containing validated binder sequences
        patient_data: Dict with patient-specific information
    
    Returns:
        Dict containing personalized results
    """
    engine = PersonalizationEngine()
    personalized_binders = []
    
    for binder in validated_input.get('validated_binders', []):
        try:
            result = engine.personalize_sequence(binder['sequence'], patient_data)
            personalized_binders.append(result)
        except Exception as e:
            logger.error(f"Failed to personalize sequence: {e}")
            continue
    
    # Sort by final score
    personalized_binders.sort(key=lambda x: x['final_score'], reverse=True)
    
    # Calculate summary statistics
    summary = {
        'total_candidates': len(personalized_binders),
        'high_celtic_scores': sum(1 for b in personalized_binders if b['celtic_score'] > 0.2),
        'average_coverage': round(
            sum(b['population_coverage'] for b in personalized_binders) / 
            len(personalized_binders) if personalized_binders else 0,
            3
        )
    }
    
    return {
        'personalized_binders': personalized_binders,
        'summary': summary,
        'metadata': {
            'ancestry_profile': patient_data.get('ancestry_profile', ['celtic_british']),
            'hla_version': engine.hla_data['metadata']['version'] if engine.hla_data else 'N/A',
            'timestamp': '2025-09-26'
        }
    }
