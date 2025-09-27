"""
Simple sequence validator for the test script.
"""
import json

class SimpleValidator:
    def __init__(self, config_path):
        with open(config_path) as f:
            self.config = json.load(f)
    
    def validate_sequence(self, sequence):
        results = {
            "valid": True,
            "warnings": [],
            "metrics": {
                "aromatic_content": 0,
                "hydrophobic_content": 0,
                "net_charge": 0,
                "binding_motifs": {"scores": {}, "total_score": 0}
            },
            "population_scores": {}
        }
        
        # Calculate metrics
        aromatics = ['F', 'W', 'Y']
        hydrophobics = ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'C']
        positives = ['R', 'K', 'H']
        negatives = ['D', 'E']
        
        seq_len = len(sequence)
        results["metrics"]["aromatic_content"] = sum(sequence.count(aa) for aa in aromatics) / seq_len * 100
        results["metrics"]["hydrophobic_content"] = sum(sequence.count(aa) for aa in hydrophobics) / seq_len * 100
        results["metrics"]["net_charge"] = (
            sum(sequence.count(aa) for aa in positives) - 
            sum(sequence.count(aa) for aa in negatives)
        )
        
        # Check against population parameters
        for pop_name, pop_params in self.config["populations"].items():
            weight = pop_params["ancestry_weight"]
            biophysical = pop_params["biophysical_params"]
            score = 1.0
            
            # Check metrics against ranges
            if not (biophysical["aromatic_content"]["min"] <= results["metrics"]["aromatic_content"] <= 
                   biophysical["aromatic_content"]["max"]):
                score *= 0.7
            
            if not (biophysical["hydrophobic_content"]["min"] <= results["metrics"]["hydrophobic_content"] <= 
                   biophysical["hydrophobic_content"]["max"]):
                score *= 0.7
            
            if not (biophysical["net_charge"]["min"] <= results["metrics"]["net_charge"] <= 
                   biophysical["net_charge"]["max"]):
                score *= 0.7
            
            # Check binding motifs
            motifs = pop_params["binding_motifs"]
            motif_matches = sum(1 for motif in motifs if motif in sequence)
            motif_score = motif_matches / len(motifs) if motifs else 0
            
            results["metrics"]["binding_motifs"]["scores"][pop_name] = {
                "score": motif_score,
                "weighted_score": motif_score * weight
            }
            
            results["population_scores"][pop_name] = {
                "score": score * (0.3 + 0.7 * motif_score),  # Binding motifs affect up to 70% of score
                "weight": weight
            }
        
        # Calculate overall weighted score
        weighted_score = sum(
            score["score"] * score["weight"]
            for score in results["population_scores"].values()
        ) / sum(score["weight"] for score in results["population_scores"].values())
        
        if weighted_score < 0.5:
            results["valid"] = False
            results["warnings"].append(f"Overall weighted score {weighted_score:.2f} below threshold")
        
        return results