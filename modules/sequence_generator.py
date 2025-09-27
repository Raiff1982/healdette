"""
Module for generating antibody sequences using ProtGPT2.
"""

from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
import random
import json

class SequenceGenerator:
    def __init__(self, config_path):
        """Initialize the sequence generator with a configuration file."""
        self.tokenizer = AutoTokenizer.from_pretrained("nferruz/ProtGPT2")
        self.model = AutoModelForCausalLM.from_pretrained("nferruz/ProtGPT2")
        self.config_path = config_path
        
        # Load Celtic-specific parameters
        with open(config_path) as f:
            config = json.load(f)
            self.celtic_params = config["populations"]["celtic"]["biophysical_params"]
            self.celtic_motifs = config["populations"]["celtic"]["binding_motifs"]

    def _check_homopolymer(self, sequence, max_repeat=4):
        """Check for homopolymer repeats."""
        for i in range(len(sequence) - max_repeat + 1):
            if len(set(sequence[i:i + max_repeat])) == 1:
                return True
        return False
    
    def _check_protein_realism(self, sequence):
        """Check protein sequence realism based on common patterns."""
        # Check for basic protein patterns
        rules = {
            'proline_patterns': len([i for i in range(len(sequence)-1) if sequence[i:i+2] == 'PP']) <= 2,
            'charged_spacing': not any(all(aa in 'KR' for aa in sequence[i:i+4]) for i in range(len(sequence)-3)),
            'hydrophobic_balance': 0.2 <= len([aa for aa in sequence if aa in 'AILMFWV'])/len(sequence) <= 0.6,
            'glycine_spacing': not any(all(aa == 'G' for aa in sequence[i:i+3]) for i in range(len(sequence)-2)),
            'cysteine_pairs': sequence.count('C') % 2 == 0,  # Cysteines should appear in pairs
            'charged_distribution': all(sequence[i:i+5].count('K') + sequence[i:i+5].count('R') <= 3 
                                     for i in range(len(sequence)-4))
        }
        return all(rules.values())

    def _get_realistic_seed(self):
        """Generate a realistic seed sequence with Celtic properties."""
        # Common antibody framework regions mixed with Celtic motifs
        framework_segments = [
            "EVQLVESGGGLVQ",  # FR1-like
            "WVRQAPGKGLEWVS",  # FR2-like
            "RFTISRDNSKNTLYLQMNSLRAEDTAVYYC"  # FR3-like
        ]
        
        # Add Celtic-specific motifs
        celtic_segments = [
            "KRHKFY",  # Charged and aromatic
            "WYRF",    # Celtic binding motif
            "KHYF"     # Celtic binding motif
        ]
        
        seed = random.choice(framework_segments) + random.choice(celtic_segments)
        return seed

    def generate_sequences(self, num_sequences=5):
        """Generate antibody sequences based on Celtic-optimized parameters with realism checks."""
        sequences = []
        max_attempts = num_sequences * 10  # Allow more attempts for quality sequences
        attempts = 0
        
        while len(sequences) < num_sequences and attempts < max_attempts:
            # Get a realistic seed sequence
            seed = self._get_realistic_seed()
            input_ids = self.tokenizer.encode(seed, return_tensors="pt")
            
            output = self.model.generate(
                input_ids,
                do_sample=True,
                top_k=50,
                top_p=0.92,
                temperature=0.8,
                max_length=120,
                min_length=60,
                pad_token_id=self.tokenizer.eos_token_id,
                repetition_penalty=1.3,
                no_repeat_ngram_size=3  # Prevent repetitive patterns
            )
            
            sequence = self.tokenizer.decode(output[0], skip_special_tokens=True)
            
            # Skip if homopolymer or unrealistic patterns found
            if self._check_homopolymer(sequence) or not self._check_protein_realism(sequence):
                attempts += 1
                continue
            
            # Calculate sequence properties
            aromatic_aas = ['F', 'W', 'Y']
            hydrophobic_aas = ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'C']
            positive_aas = ['R', 'K', 'H']
            negative_aas = ['D', 'E']
            
            seq_len = len(sequence)
            aromatic_content = sum(sequence.count(aa) for aa in aromatic_aas) / seq_len * 100
            hydrophobic_content = sum(sequence.count(aa) for aa in hydrophobic_aas) / seq_len * 100
            net_charge = (sum(sequence.count(aa) for aa in positive_aas) - 
                         sum(sequence.count(aa) for aa in negative_aas))
            
            # Check Celtic criteria with protein realism
            criteria_met = 0
            if (self.celtic_params["aromatic_content"]["min"] <= aromatic_content <= 
                self.celtic_params["aromatic_content"]["max"]):
                criteria_met += 1
            if (self.celtic_params["hydrophobic_content"]["min"] <= hydrophobic_content <= 
                self.celtic_params["hydrophobic_content"]["max"]):
                criteria_met += 1
            if (self.celtic_params["net_charge"]["min"] <= net_charge <= 
                self.celtic_params["net_charge"]["max"]):
                criteria_met += 1
            
            # Accept sequences meeting minimum criteria
            meets_realism = self._check_protein_realism(sequence)
            no_homopolymer = not self._check_homopolymer(sequence)
            
            quality_score = criteria_met + (1 if meets_realism else 0) + (1 if no_homopolymer else 0)
            
            # Accept if good quality or running out of attempts
            if quality_score >= 3 or (attempts > max_attempts * 0.8 and not sequences):
                sequences.append(sequence)
            
            attempts += 1
            
            # Break early if we have enough sequences
            if len(sequences) >= num_sequences:
                break
            
            attempts += 1
            
            attempts += 1
        
        return sequences