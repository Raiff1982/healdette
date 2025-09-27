
from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
import random

"""
Module for validating antibody sequences against multi-ethnic parameters.
"""
import json
from pathlib import Path
import re

class WeightedSequenceValidator:
    def __init__(self, config_path):
        """Initialize validator with configuration file."""
        self.config_path = config_path
        self.config = self._load_config()
        
    def _load_config(self):
        """Load and parse configuration file."""
        with open(self.config_path) as f:
            return json.load(f)
    
    def _calculate_aromatic_content(self, sequence):
        """Calculate percentage of aromatic amino acids."""
        aromatics = ['F', 'W', 'Y']
        aromatic_count = sum(sequence.count(aa) for aa in aromatics)
        return (aromatic_count / len(sequence)) * 100
    
    def _calculate_hydrophobic_content(self, sequence):
        """Calculate percentage of hydrophobic amino acids."""
        hydrophobic = ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'C']
        hydrophobic_count = sum(sequence.count(aa) for aa in hydrophobic)
        return (hydrophobic_count / len(sequence)) * 100
    
    def _calculate_net_charge(self, sequence):
        """Calculate net charge of the sequence."""
        positive = ['R', 'K', 'H']
        negative = ['D', 'E']
        pos_count = sum(sequence.count(aa) for aa in positive)
        neg_count = sum(sequence.count(aa) for aa in negative)
        return pos_count - neg_count
    
    def _check_binding_motifs(self, sequence, population_motifs):
        """Check for presence of population-specific binding motifs."""
        found_motifs = 0
        total_motifs = len(population_motifs)
        
        for motif in population_motifs:
            if motif in sequence:
                found_motifs += 1
                
        return found_motifs / total_motifs if total_motifs > 0 else 0
    
    def validate_sequence(self, sequence):
        """Validate a sequence against all population parameters."""
        results = {
            "valid": True,
            "warnings": [],
            "metrics": {
                "aromatic_content": self._calculate_aromatic_content(sequence),
                "hydrophobic_content": self._calculate_hydrophobic_content(sequence),
                "net_charge": self._calculate_net_charge(sequence),
                "binding_motifs": {
                    "scores": {},
                    "total_score": 0
                }
            },
            "population_scores": {}
        }
        
        # Check sequence length
        seq_len = len(sequence)
        if not (self.config["global_params"]["sequence_length"]["min"] <= seq_len <= 
                self.config["global_params"]["sequence_length"]["max"]):
            results["warnings"].append(f"Sequence length {seq_len} outside acceptable range")
            results["valid"] = False
        
        # Validate against each population's parameters
        total_weight = 0
        weighted_score = 0
        
        for pop_name, pop_params in self.config["populations"].items():
            pop_score = 1.0  # Start with perfect score
            weight = pop_params["ancestry_weight"]
            biophysical = pop_params["biophysical_params"]
            
            # Check aromatic content
            if not (biophysical["aromatic_content"]["min"] <= results["metrics"]["aromatic_content"] <= 
                   biophysical["aromatic_content"]["max"]):
                pop_score *= 0.7
            
            # Check hydrophobic content
            if not (biophysical["hydrophobic_content"]["min"] <= results["metrics"]["hydrophobic_content"] <= 
                   biophysical["hydrophobic_content"]["max"]):
                pop_score *= 0.7
            
            # Check net charge
            if not (biophysical["net_charge"]["min"] <= results["metrics"]["net_charge"] <= 
                   biophysical["net_charge"]["max"]):
                pop_score *= 0.7
            
            # Check binding motifs
            motif_score = self._check_binding_motifs(sequence, pop_params["binding_motifs"])
            results["metrics"]["binding_motifs"]["scores"][pop_name] = {
                "score": motif_score,
                "weighted_score": motif_score * weight
            }
            pop_score *= (0.3 + 0.7 * motif_score)  # Binding motifs affect up to 70% of score
            
            results["population_scores"][pop_name] = {
                "score": pop_score,
                "weight": weight
            }
            
            total_weight += weight
            weighted_score += pop_score * weight
        
        # Calculate total binding motif score
        results["metrics"]["binding_motifs"]["total_score"] = sum(
            score["weighted_score"] 
            for score in results["metrics"]["binding_motifs"]["scores"].values()
        )
        
        # Normalize weighted score
        if total_weight > 0:
            weighted_score /= total_weight
        
        # Update validity based on overall score
        if weighted_score < 0.5:
            results["valid"] = False
            results["warnings"].append(f"Overall weighted score {weighted_score:.2f} below threshold")
        
        return results
        max_length=200,
        num_return_sequences=num_candidates
    )

    binders = []
    for output in outputs:
        sequence = tokenizer.decode(output, skip_special_tokens=True)
        sequence = ''.join([aa for aa in sequence if aa in "ACDEFGHIKLMNPQRSTVWY"])
        if len(sequence) > 30:
            binder_meta = {
                "sequence": sequence,
                "perspective_source": fusion_context["perspective_tags"],
                "sentiment_trace": fusion_context["sentiment_trace"],
                "symbolic_logic_score": fusion_context["symbolic_logic_score"]
            }
            binders.append(binder_meta)

    return {"generated_binders": binders}
