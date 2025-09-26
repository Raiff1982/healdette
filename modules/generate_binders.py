
from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
import random
from typing import List, Dict
import json
import os

class AntibodyGenerator:
    # Typical antibody framework regions for sequence templating
    FRAMEWORK_REGIONS = {
        'VH': {
            'FR1': 'EVQLVESGGGLVQPGGSLRLSCAAS',
            'FR2': 'WVRQAPGKGLEWVS',
            'FR3': 'RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR',
            'FR4': 'WGQGTLVTVSS'
        },
        'VL': {
            'FR1': 'DIQMTQSPSSLSASVGDRVTITC',
            'FR2': 'WYQQKPGKAPKLLIY',
            'FR3': 'GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC',
            'FR4': 'FGQGTKVEIK'
        }
    }

    def __init__(self):
        """Initialize the antibody generator with ProtGPT2 model and validation data."""
        # Load ProtGPT2 model
        self.tokenizer = AutoTokenizer.from_pretrained("nferruz/ProtGPT2")
        self.model = AutoModelForCausalLM.from_pretrained("nferruz/ProtGPT2")
        
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
            self.model.config.pad_token_id = self.model.config.eos_token_id

        # Load validation dataset of known antibodies
        self.validation_set = self._load_validation_set()

    def _load_validation_set(self) -> List[Dict]:
        """Load curated set of known therapeutic antibodies for validation."""
        validation_file = os.path.join(os.path.dirname(__file__), 
                                     'data', 'therapeutic_antibodies.json')
        try:
            with open(validation_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print("Warning: Validation dataset not found. Using basic validation.")
            return []

    def _generate_cdrs(self, context: str, num_variants: int = 5) -> List[str]:
        """Generate CDR sequences with antibody-specific constraints."""
        inputs = self.tokenizer(context, return_tensors="pt", padding=True)
        
        outputs = self.model.generate(
            inputs["input_ids"],
            attention_mask=inputs["attention_mask"],
            do_sample=True,
            top_k=50,  # More focused sampling
            top_p=0.95,
            temperature=0.8,  # Slightly lower temperature for more conservative sampling
            max_length=50,
            num_return_sequences=num_variants,
            pad_token_id=self.tokenizer.pad_token_id,
            no_repeat_ngram_size=3  # Prevent repetitive motifs
        )
        
        cdrs = []
        for output in outputs:
            sequence = self.tokenizer.decode(output, skip_special_tokens=True)
            # Clean and validate CDR sequence
            sequence = ''.join(aa for aa in sequence if aa in "ACDEFGHIKLMNPQRSTVWY")
            if 5 <= len(sequence) <= 20:  # Typical CDR length constraints
                cdrs.append(sequence)
        
        return cdrs

    def _assemble_antibody(self, heavy_cdrs: List[str], light_cdrs: List[str]) -> str:
        """Assemble complete antibody sequence from framework regions and CDRs."""
        vh = (self.FRAMEWORK_REGIONS['VH']['FR1'] + 
              heavy_cdrs[0] + 
              self.FRAMEWORK_REGIONS['VH']['FR2'] + 
              heavy_cdrs[1] + 
              self.FRAMEWORK_REGIONS['VH']['FR3'] + 
              heavy_cdrs[2] + 
              self.FRAMEWORK_REGIONS['VH']['FR4'])
        
        vl = (self.FRAMEWORK_REGIONS['VL']['FR1'] + 
              light_cdrs[0] + 
              self.FRAMEWORK_REGIONS['VL']['FR2'] + 
              light_cdrs[1] + 
              self.FRAMEWORK_REGIONS['VL']['FR3'] + 
              light_cdrs[2] + 
              self.FRAMEWORK_REGIONS['VL']['FR4'])
        
        return vh + vl

    def generate_binders(self, fusion_context: Dict, num_candidates: int = 10) -> Dict:
        """Generate antibody sequences using template-based approach with ProtGPT2."""
        binders = []
        
        # Extract targeting information from fusion context
        target_motif = fusion_context.get('cleaned_sequence', '')[:10]
        
        for _ in range(num_candidates):
            # Generate CDRs with target context
            heavy_cdrs = self._generate_cdrs(f"Target binding site: {target_motif}", 3)
            light_cdrs = self._generate_cdrs(f"Light chain CDRs for {target_motif}", 3)
            
            if len(heavy_cdrs) >= 3 and len(light_cdrs) >= 3:
                sequence = self._assemble_antibody(heavy_cdrs, light_cdrs)
                binders.append({
                    "sequence": sequence,
                    "heavy_cdrs": heavy_cdrs[:3],
                    "light_cdrs": light_cdrs[:3],
                    "target_context": target_motif
                })
        
        return {"generated_binders": binders}

    return {"generated_binders": binders}
