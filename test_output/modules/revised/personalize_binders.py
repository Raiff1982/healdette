"""Personalize binders based on patient HLA types and ancestry."""

from typing import Dict, List, Optional
from pathlib import Path
import json

# HLA reference profiles based on ancestry
HLA_REFERENCE = {
    "Celtic_British": ["A*01:01", "A*03:01", "B*07:02", "B*08:01", "C*07:01"],
    "Irish_Dominant": ["A*01:01", "B*08:01", "C*07:01"],
    "Scottish": ["A*02:01", "B*07:02", "C*04:01"],
    "Germanic": ["A*02:01", "A*03:01", "B*07:02", "C*07:01"]
}

class BiomarkerProfile:
    """Patient biomarker and ancestry profile."""
    
    def __init__(self, hla_types: List[str], ancestry_weights: Dict[str, float]):
        """Initialize biomarker profile.
        
        Args:
            hla_types: List of patient HLA types (e.g., ["A*01:01", "B*08:01"])
            ancestry_weights: Dict mapping ancestry to weight (e.g., {"irish": 0.42})
        """
        self.hla_types = hla_types
        self.ancestry_weights = ancestry_weights
        self.hla_frequencies = self._load_hla_frequencies()
        
    def _load_hla_frequencies(self) -> Dict:
        """Load population HLA frequency data."""
        data_file = Path(__file__).parent / "data" / "hla_frequencies.json"
        with open(data_file) as f:
            return json.load(f)["populations"]
            
    def get_hla_compatibility_score(self, sequence: str) -> float:
        """Calculate HLA compatibility score for a sequence.
        
        Args:
            sequence: Amino acid sequence to evaluate
            
        Returns:
            float: Compatibility score between 0-1
        """
        # Calculate weighted score based on patient HLA types and ancestry
        # This is a simplified scoring function - in practice would use more
        # sophisticated epitope prediction
        
        celtic_markers = set(["A*03:01", "B*08:01", "C*07:01"])
        british_markers = set(["A*02:01", "B*07:02", "C*04:01"])
        
        score = 0.0
        n_matches = 0
        
        # Weight by ancestry profile
        for hla in self.hla_types:
            if hla in celtic_markers:
                score += (self.ancestry_weights.get("irish", 0) * 0.4 + 
                         self.ancestry_weights.get("scottish", 0) * 0.3)
                n_matches += 1
            if hla in british_markers:
                score += (self.ancestry_weights.get("english", 0) * 0.4 +
                         self.ancestry_weights.get("welsh", 0) * 0.2)
                n_matches += 1
                
        # Normalize score
        return score / max(1, n_matches)
    
    def personalize_binder(self, binder: Dict) -> Dict:
        """Personalize a binder based on patient profile.
        
        Args:
            binder: Dict containing binder sequence and properties
            
        Returns:
            Dict: Binder with added personalization scores
        """
        sequence = binder["sequence"]
        
        # Calculate HLA compatibility
        hla_score = self.get_hla_compatibility_score(sequence)
        
        # Add personalization scores
        binder["personalization"] = {
            "hla_compatibility": hla_score,
            "matched_hla_types": [hla for hla in self.hla_types 
                                if hla in HLA_REFERENCE["Celtic_British"]],
            "ancestry_profile": {
                "primary": "Celtic_British",
                "weights": self.ancestry_weights
            }
        }
        
        return binder
        
def personalize_binders(binders: List[Dict], patient_data: Dict) -> List[Dict]:
    """Personalize a list of binders for a specific patient.
    
    Args:
        binders: List of binder dicts with sequences and properties
        patient_data: Dict containing patient HLA types and ancestry info
        
    Returns:
        List[Dict]: Personalized binders with compatibility scores
    """
    profile = BiomarkerProfile(
        hla_types=patient_data["immune_profile"],
        ancestry_weights=patient_data["population_weights"]
    )
    
    personalized = []
    for binder in binders:
        personalized.append(profile.personalize_binder(binder))
        
    # Sort by compatibility score
    personalized.sort(key=lambda x: x["personalization"]["hla_compatibility"], 
                     reverse=True)
        
    return personalized