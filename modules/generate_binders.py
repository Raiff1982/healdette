
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
        """Initialize the antibody sequence generator.
        
        Sets up the ProtGPT2 language model for sequence generation and loads the validation dataset.
        The model is configured with padding tokens and loaded with pre-trained weights from
        nferruz/ProtGPT2, which was trained on antibody sequences.
        """
        # Initialize ProtGPT2 model and tokenizer
        self.tokenizer = AutoTokenizer.from_pretrained("nferruz/ProtGPT2")
        self.model = AutoModelForCausalLM.from_pretrained("nferruz/ProtGPT2")
        
        # Configure padding tokens for batch processing
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
            self.model.config.pad_token_id = self.model.config.eos_token_id

        # Load validation dataset from THERAb database
        self.validation_set = self._load_validation_set()

    def _load_validation_set(self) -> List[Dict]:
        """Load the validation dataset of therapeutic antibodies.
        
        Returns:
            List[Dict]: A list of validated therapeutic antibodies, each containing:
                - sequence: Full amino acid sequence
                - name: Antibody name/identifier
                - target: Target antigen
                - approved: FDA approval status
        """
        validation_file = os.path.join(os.path.dirname(__file__), 
                                     'data', 'therapeutic_antibodies.json')
        try:
            with open(validation_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print("Warning: Validation dataset not found. Using basic validation.")
            return []

    def _generate_cdrs(self, context: str, num_variants: int = 5) -> List[str]:
        """Generate complementarity determining region (CDR) sequences.
        
        Uses ProtGPT2 to generate CDR sequences with antibody-specific constraints:
        - Valid amino acid composition
        - Length constraints (5-20 residues)
        - No repetitive motifs (using n-gram blocking)
        
        Args:
            context (str): Input sequence context for conditional generation
            num_variants (int, optional): Number of sequences to generate. Defaults to 5.
        
        Returns:
            List[str]: List of valid CDR sequences meeting constraints
        """
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
        """Assemble a complete antibody sequence by combining framework regions and CDRs.
        
        Constructs the full antibody sequence by inserting CDR sequences into their
        corresponding positions within the conserved framework regions. Uses standard
        IMGT-defined framework regions for VH and VL domains.
        
        Args:
            heavy_cdrs (List[str]): Three heavy chain CDR sequences [CDR-H1, CDR-H2, CDR-H3]
            light_cdrs (List[str]): Three light chain CDR sequences [CDR-L1, CDR-L2, CDR-L3]
            
        Returns:
            str: Complete antibody sequence with framework regions and CDRs
        """
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
        """Generate a set of candidate antibody sequences for a given target.
        
        Uses ProtGPT2 to generate CDR sequences conditioned on the target sequence,
        then assembles complete antibodies using template-based framework regions.
        Each candidate is validated against known therapeutic antibodies.
        
        Args:
            fusion_context (Dict): Context for sequence generation containing:
                - cleaned_sequence (str): Target protein sequence
                - optional additional parameters
            num_candidates (int, optional): Number of sequences to generate. Defaults to 10.
            
        Returns:
            Dict: Generated antibody sequences and metadata
                - generated_binders (List[Dict]): List of candidate antibodies
                    - sequence (str): Full antibody sequence
                    - heavy_cdrs (List[str]): Heavy chain CDR sequences
                    - light_cdrs (List[str]): Light chain CDR sequences
                    - validation_score (float): Similarity to known antibodies
        """
        binders = []
        
        # Extract targeting information from fusion context
        target_motif = fusion_context.get('cleaned_sequence', '')[:10]
        
        for _ in range(num_candidates):
            # Generate CDRs with target context
            heavy_cdrs = self._generate_cdrs(f"Target binding site: {target_motif}", 3)
            light_cdrs = self._generate_cdrs(f"Light chain CDRs for {target_motif}", 3)
            
            if len(heavy_cdrs) >= 3 and len(light_cdrs) >= 3:
                sequence = self._assemble_antibody(heavy_cdrs[:3], light_cdrs[:3])
                validation_score = self._validate_sequence(sequence)
                binder = {
                    "sequence": sequence,
                    "heavy_cdrs": heavy_cdrs[:3],
                    "light_cdrs": light_cdrs[:3],
                    "validation_score": validation_score
                }
                binders.append(binder)
                    "sequence": sequence,
                    "heavy_cdrs": heavy_cdrs,
                    "light_cdrs": light_cdrs,
                    "validation_score": self._validate_sequence(sequence)
                })
        
        return {"generated_binders": binders}
    
    def _validate_sequence(self, sequence: str) -> float:
        """Validate a generated sequence against known therapeutic antibodies.
        
        Compares the generated sequence with validated therapeutic antibodies
        to assess its similarity to known successful designs.
        
        Args:
            sequence (str): Antibody sequence to validate
            
        Returns:
            float: Validation score between 0 and 1
        """
        if not self.validation_set:
            return 0.5  # Default score when validation set is unavailable
            
        max_score = 0
        for ref_ab in self.validation_set:
            score = self._calculate_similarity(sequence, ref_ab["sequence"])
            max_score = max(max_score, score)
            
        return max_score
                    "sequence": sequence,
                    "heavy_cdrs": heavy_cdrs[:3],
                    "light_cdrs": light_cdrs[:3],
                    "target_context": target_motif
                })
        
        return {"generated_binders": binders}

    return {"generated_binders": binders}
