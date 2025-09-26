"""Antibody sequence generation using ProtGPT2 with IMGT germline templates."""

from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
import random
from typing import List, Dict
import json
import os

class AntibodyGenerator:
    """Generates antibody sequences using ProtGPT2 with IMGT germline templates."""
    
    # IMGT human germline V-region templates
    GERMLINE_TEMPLATES = {
        'VH': {
            'IGHV1-69*01': {
                'FR1': 'QVQLVQSGAEVKKPGSSVKVSCKASGGTFS',
                'FR2': 'WVRQAPGQGLEWMG',
                'FR3': 'RVTITADKSTSTAYMELSSLRSEDTAVYYCAR',
                'FR4': 'WGQGTLVTVSS'
            },
            'IGHV3-23*01': {
                'FR1': 'EVQLLESGGGLVQPGGSLRLSCAASGFTFS',
                'FR2': 'WVRQAPGKGLEWVS',
                'FR3': 'RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR',
                'FR4': 'WGQGTLVTVSS'
            }
        },
        'VL': {
            'IGKV1-39*01': {
                'FR1': 'DIQMTQSPSSLSASVGDRVTITC',
                'FR2': 'WYQQKPGKAPKLLIY',
                'FR3': 'GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC',
                'FR4': 'FGQGTKVEIK'
            },
            'IGLV1-44*01': {
                'FR1': 'QSVLTQPPSVSGAPGQRVTISCTGSSSNIG',
                'FR2': 'WYQQHPGKAPKLLIY',
                'FR3': 'GVPDRFSGSGSGTDFTLTISGVQAEDVAVYYC',
                'FR4': 'FGGGTKLTVL'
            }
        }
    }
    
    # Sequence quality parameters
    MAX_HOMOPOLYMER_LENGTH = 4
    MIN_CDR_LENGTH = 5
    MAX_CDR_LENGTH = 20
    MIN_DISORDER_SCORE = 0.0
    MAX_DISORDER_SCORE = 0.3
    ALLOWED_PI_RANGE = (5.5, 8.5)

    def __init__(self):
        """Initialize the antibody sequence generator."""
        self.tokenizer = AutoTokenizer.from_pretrained("nferruz/ProtGPT2")
        self.model = AutoModelForCausalLM.from_pretrained("nferruz/ProtGPT2")
        
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
            self.model.config.pad_token_id = self.model.config.eos_token_id

        self.validation_set = self._load_validation_set()

    def _load_validation_set(self) -> List[Dict]:
        """Load the validation dataset of therapeutic antibodies."""
        validation_file = os.path.join(os.path.dirname(__file__), 
                                     'data', 'therapeutic_antibodies.json')
        try:
            with open(validation_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print("Warning: Validation dataset not found. Using basic validation.")
            return []

    def _generate_cdrs(self, context: str, template_seq: str, num_variants: int = 5) -> List[str]:
        """Generate CDR sequences with template-based conditioning."""
        prompt = f"{template_seq} {context} <CDR>"
        inputs = self.tokenizer(prompt, return_tensors="pt", padding=True)
        
        outputs = self.model.generate(
            inputs["input_ids"],
            attention_mask=inputs["attention_mask"],
            do_sample=True,
            top_k=40,  # Allow more diversity
            top_p=0.9,  # Slightly more permissive
            temperature=0.7,  # Slightly higher for more diversity
            max_new_tokens=25,  # Allow slightly longer sequences
            num_return_sequences=num_variants * 3,  # Generate more candidates
            pad_token_id=self.tokenizer.pad_token_id,
            no_repeat_ngram_size=2,
            repetition_penalty=1.3  # Slightly more permissive
        )
        
        cdrs = []
        for output in outputs:
            sequence = self.tokenizer.decode(output, skip_special_tokens=True)
            sequence = sequence.split("<CDR>")[-1].strip()
            sequence = ''.join(aa for aa in sequence if aa in "ACDEFGHIKLMNPQRSTVWY")
            
            if self._check_sequence_quality(sequence):
                cdrs.append(sequence)
            
            if len(cdrs) >= num_variants:
                break
                
        return cdrs

    def _check_sequence_quality(self, sequence: str) -> bool:
        """Check if a sequence meets all quality criteria."""
        if not (self.MIN_CDR_LENGTH <= len(sequence) <= self.MAX_CDR_LENGTH):
            return False
            
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            if aa * self.MAX_HOMOPOLYMER_LENGTH in sequence:
                return False
                
        aa_counts = {aa: sequence.count(aa) / len(sequence) for aa in set(sequence)}
        if max(aa_counts.values()) > 0.3:
            return False
            
        hydrophobic = "AILMFWYV"
        hydrophobic_fraction = sum(sequence.count(aa) for aa in hydrophobic) / len(sequence)
        if hydrophobic_fraction > 0.6 or hydrophobic_fraction < 0.2:
            return False
            
        return True

    def _assemble_antibody(self, heavy_cdrs: List[str], light_cdrs: List[str]) -> str:
        """Assemble complete antibody sequence from framework regions and CDRs."""
        vh = (self.GERMLINE_TEMPLATES['VH']['IGHV3-23*01']['FR1'] + 
              heavy_cdrs[0] + 
              self.GERMLINE_TEMPLATES['VH']['IGHV3-23*01']['FR2'] + 
              heavy_cdrs[1] + 
              self.GERMLINE_TEMPLATES['VH']['IGHV3-23*01']['FR3'] + 
              heavy_cdrs[2] + 
              self.GERMLINE_TEMPLATES['VH']['IGHV3-23*01']['FR4'])
        
        vl = (self.GERMLINE_TEMPLATES['VL']['IGKV1-39*01']['FR1'] + 
              light_cdrs[0] + 
              self.GERMLINE_TEMPLATES['VL']['IGKV1-39*01']['FR2'] + 
              light_cdrs[1] + 
              self.GERMLINE_TEMPLATES['VL']['IGKV1-39*01']['FR3'] + 
              light_cdrs[2] + 
              self.GERMLINE_TEMPLATES['VL']['IGKV1-39*01']['FR4'])
        
        return vh + vl

    def generate_binders(self, fusion_context: Dict, num_candidates: int = 10) -> Dict:
        """Generate a set of candidate antibody sequences for a given target."""
        binders = []
        target_motif = fusion_context.get('cleaned_sequence', '')[:20]
        
        vh_templates = list(self.GERMLINE_TEMPLATES['VH'].items())
        vl_templates = list(self.GERMLINE_TEMPLATES['VL'].items())
        
        attempts = 0
        max_attempts = num_candidates * 3
        
        while len(binders) < num_candidates and attempts < max_attempts:
            attempts += 1
            
            vh_name, vh = random.choice(vh_templates)
            vl_name, vl = random.choice(vl_templates)
            
            heavy_cdrs = self._generate_cdrs(
                f"Target binding site: {target_motif}",
                vh['FR1'] + "X" * 10 + vh['FR2'],
                3
            )
            
            light_cdrs = self._generate_cdrs(
                f"Light chain CDRs for {target_motif}",
                vl['FR1'] + "X" * 8 + vl['FR2'],
                3
            )
            
            if len(heavy_cdrs) >= 3 and len(light_cdrs) >= 3:
                sequence = self._assemble_antibody(heavy_cdrs[:3], light_cdrs[:3])
                validation_score = self._validate_sequence(sequence)
                
                if validation_score >= 0.7:
                    binder = {
                        "sequence": sequence,
                        "heavy_cdrs": heavy_cdrs[:3],
                        "light_cdrs": light_cdrs[:3],
                        "validation_score": validation_score,
                        "template_vh": vh_name,
                        "template_vl": vl_name
                    }
                    binders.append(binder)
        
        return {
            "generated_binders": binders,
            "stats": {
                "attempts": attempts,
                "success_rate": len(binders) / attempts if attempts > 0 else 0
            }
        }

    def _validate_sequence(self, sequence: str) -> float:
        """Validate a sequence against known therapeutic antibodies."""
        if not self.validation_set:
            return 0.5
            
        max_score = 0
        for ref_ab in self.validation_set:
            score = self._calculate_similarity(sequence, ref_ab["sequence"])
            max_score = max(max_score, score)
            
        return max_score

    def _calculate_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity between two antibody sequences."""
        min_len = min(len(seq1), len(seq2))
        matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
        return matches / min_len if min_len > 0 else 0.0