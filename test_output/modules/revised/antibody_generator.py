"""Antibody sequence generation using ProtGPT2 with IMGT germline templates."""

from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
import random
from typing import List, Dict, Tuple
import json
import os
from pathlib import Path
from sequence_validator import SequenceValidator

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

    def __init__(self):
        """Initialize generator with model and validator."""
        # Initialize ProtGPT2
        self.tokenizer = AutoTokenizer.from_pretrained("nferruz/ProtGPT2")
        self.model = AutoModelForCausalLM.from_pretrained("nferruz/ProtGPT2")
        
        # Configure tokenizer
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
            self.model.config.pad_token_id = self.model.config.eos_token_id
        
        # Initialize sequence validator
        self.validator = SequenceValidator()
        
        # Load validation dataset
        self.validation_set = self._load_validation_set()

    def _load_validation_set(self) -> List[Dict]:
        """Load therapeutic antibody validation dataset."""
        try:
            validation_file = Path(__file__).parent / 'data' / 'therapeutic_antibodies.json'
            with open(validation_file) as f:
                data = json.load(f)
                return data.get('therapeutic_antibodies', [])
        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Warning: Error loading validation dataset: {str(e)}. Using basic validation.")
            return []

    def _generate_cdrs(self, context: str, template_seq: str, num_variants: int = 5) -> List[Tuple[str, Dict]]:
        """Generate CDR sequences with template-based conditioning."""
        prompt = f"{template_seq} {context} <CDR>"
        inputs = self.tokenizer(prompt, return_tensors="pt", padding=True)
        
        print(f"Generating CDRs with prompt: {prompt}")
        
        # Generate sequences with improved parameters
        outputs = self.model.generate(
            inputs["input_ids"],
            attention_mask=inputs["attention_mask"],
            do_sample=True,
            top_k=50,
            top_p=0.95,
            temperature=0.8,
            max_new_tokens=30,
            min_new_tokens=15,
            num_return_sequences=num_variants * 5,
            pad_token_id=self.tokenizer.pad_token_id,
            no_repeat_ngram_size=3,
            repetition_penalty=1.5,
            bad_words_ids=[[self.tokenizer.eos_token_id]]
        )
        
        # Process and validate sequences
        valid_sequences = []
        for i, output in enumerate(outputs):
            sequence = self.tokenizer.decode(output, skip_special_tokens=True)
            sequence = sequence.split("<CDR>")[-1].strip()
            sequence = ''.join(aa for aa in sequence if aa in "ACDEFGHIKLMNPQRSTVWY")
            
            print(f"Generated sequence {i+1}: {sequence}")
            
            # Analyze sequence quality
            analysis = self.validator.analyze_sequence(sequence)
            print(f"Validation results for sequence {i+1}:")
            print(f"  Passes validation: {analysis['passes_validation']}")
            print(f"  Analysis: {analysis}")
            
            if analysis['passes_validation']:
                valid_sequences.append((sequence, analysis))
                if len(valid_sequences) >= num_variants:
                    break
        
        return valid_sequences

    def _assemble_antibody(self, heavy_cdrs: List[str], light_cdrs: List[str], vh_template: Dict, vl_template: Dict) -> str:
        """Assemble complete antibody sequence from CDRs and templates."""
        vh = (vh_template['FR1'] + 
              heavy_cdrs[0] + 
              vh_template['FR2'] + 
              heavy_cdrs[1] + 
              vh_template['FR3'] + 
              heavy_cdrs[2] + 
              vh_template['FR4'])
        
        vl = (vl_template['FR1'] + 
              light_cdrs[0] + 
              vl_template['FR2'] + 
              light_cdrs[1] + 
              vl_template['FR3'] + 
              light_cdrs[2] + 
              vl_template['FR4'])
        
        return vh + vl

    def generate_binders(self, fusion_context: Dict, num_candidates: int = 10) -> Dict:
        """Generate antibody sequences with comprehensive validation."""
        binders = []
        target_motif = fusion_context.get('cleaned_sequence', '')[:20]
        
        vh_templates = list(self.GERMLINE_TEMPLATES['VH'].items())
        vl_templates = list(self.GERMLINE_TEMPLATES['VL'].items())
        
        attempts = 0
        max_attempts = num_candidates * 3
        
        while len(binders) < num_candidates and attempts < max_attempts:
            attempts += 1
            
            # Select templates
            vh_name, vh = random.choice(vh_templates)
            vl_name, vl = random.choice(vl_templates)
            
            # Generate and validate CDRs
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
                # Use only the sequences, not their analysis dicts
                h_seqs = [seq for seq, _ in heavy_cdrs[:3]]
                l_seqs = [seq for seq, _ in light_cdrs[:3]]
                
                sequence = self._assemble_antibody(h_seqs, l_seqs, vh, vl)
                validation_score = self._validate_sequence(sequence)
                
                if validation_score >= 0.5:
                    binder = {
                        "sequence": sequence,
                        "heavy_cdrs": h_seqs,
                        "light_cdrs": l_seqs,
                        "validation_score": validation_score,
                        "template_vh": vh_name,
                        "template_vl": vl_name,
                        "heavy_cdr_analysis": [analysis for _, analysis in heavy_cdrs[:3]],
                        "light_cdr_analysis": [analysis for _, analysis in light_cdrs[:3]]
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
        """Validate sequence against known therapeutic antibodies."""
        if not self.validation_set:
            return 0.5
        
        max_score = 0
        for ref_ab in self.validation_set:
            score = self._calculate_similarity(sequence, ref_ab["sequence"])
            max_score = max(max_score, score)
        
        return max_score

    def _calculate_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity using local alignment."""
        min_len = min(len(seq1), len(seq2))
        matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
        return matches / min_len if min_len > 0 else 0.0