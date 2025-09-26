
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
    MAX_HOMOPOLYMER_LENGTH = 4  # Maximum allowed length of amino acid repeats
    MIN_CDR_LENGTH = 5  # Minimum CDR length
    MAX_CDR_LENGTH = 20  # Maximum CDR length
    MIN_DISORDER_SCORE = 0.0  # Minimum allowed disorder score
    MAX_DISORDER_SCORE = 0.3  # Maximum allowed disorder score
    ALLOWED_PI_RANGE = (5.5, 8.5)  # Allowed range for isoelectric point

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

    def _generate_cdrs(self, context: str, template_seq: str, num_variants: int = 5) -> List[str]:
        """Generate complementarity determining region (CDR) sequences.
        
        Uses ProtGPT2 to generate CDR sequences with template-based conditioning and
        strict quality controls to prevent unrealistic sequences.
        
        Args:
            context (str): Target binding context for generation
            template_seq (str): Template sequence from germline to condition on
            num_variants (int, optional): Number of sequences to generate. Defaults to 5.
        
        Returns:
            List[str]: List of valid CDR sequences meeting all quality criteria
        """
        # Prepare input with both context and template
        prompt = f"{template_seq} {context} <CDR>"
        inputs = self.tokenizer(prompt, return_tensors="pt", padding=True)
        
        # Use more conservative sampling parameters
        outputs = self.model.generate(
            inputs["input_ids"],
            attention_mask=inputs["attention_mask"],
            do_sample=True,
            top_k=20,  # More restrictive top-k
            top_p=0.85,  # More conservative nucleus sampling
            temperature=0.6,  # Lower temperature for more conservative sampling
            max_length=30,  # Shorter max length for CDRs
            num_return_sequences=num_variants * 2,  # Generate extra for filtering
            pad_token_id=self.tokenizer.pad_token_id,
            no_repeat_ngram_size=2,  # Stricter repeat prevention
            repetition_penalty=1.5,  # Additional penalty for repetition
            length_penalty=0.8  # Slight penalty for longer sequences
        )
        
        cdrs = []
        for output in outputs:
            sequence = self.tokenizer.decode(output, skip_special_tokens=True)
            sequence = sequence.split("<CDR>")[-1].strip()  # Extract CDR part
            sequence = ''.join(aa for aa in sequence if aa in "ACDEFGHIKLMNPQRSTVWY")
            
            # Apply quality filters
            if self._check_sequence_quality(sequence):
                cdrs.append(sequence)
            
            if len(cdrs) >= num_variants:
                break
                
        return cdrs
        
    def _check_sequence_quality(self, sequence: str) -> bool:
        """Check if a sequence meets all quality criteria.
        
        Validates sequence length, amino acid composition, and structural properties.
        
        Args:
            sequence (str): Amino acid sequence to check
            
        Returns:
            bool: True if sequence meets all criteria, False otherwise
        """
        # Length check
        if not (self.MIN_CDR_LENGTH <= len(sequence) <= self.MAX_CDR_LENGTH):
            return False
            
        # Check for homopolymer runs
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            if aa * self.MAX_HOMOPOLYMER_LENGTH in sequence:
                return False
                
        # Basic composition checks
        aa_counts = {aa: sequence.count(aa) / len(sequence) for aa in set(sequence)}
        if max(aa_counts.values()) > 0.3:  # No amino acid should be >30% of sequence
            return False
            
        # Hydrophobicity balance
        hydrophobic = "AILMFWYV"
        hydrophobic_fraction = sum(sequence.count(aa) for aa in hydrophobic) / len(sequence)
        if hydrophobic_fraction > 0.6 or hydrophobic_fraction < 0.2:
            return False
            
        return True
        
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
                    - sequence (str): Full antibody sequence
                    - heavy_cdrs (List[str]): Heavy chain CDR sequences
                    - light_cdrs (List[str]): Light chain CDR sequences
                    - validation_score (float): Similarity to known antibodies
                    - template_vh (str): Heavy chain template used
                    - template_vl (str): Light chain template used
        """
        binders = []
        target_motif = fusion_context.get('cleaned_sequence', '')[:20]  # Use longer context
        
        # Select templates
        vh_templates = list(self.GERMLINE_TEMPLATES['VH'].items())
        vl_templates = list(self.GERMLINE_TEMPLATES['VL'].items())
        
        attempts = 0
        max_attempts = num_candidates * 3  # Allow multiple attempts per candidate
        
        while len(binders) < num_candidates and attempts < max_attempts:
            attempts += 1
            
            # Randomly select templates
            vh_name, vh = random.choice(vh_templates)
            vl_name, vl = random.choice(vl_templates)
            
            # Generate CDRs using templates for conditioning
            heavy_cdrs = self._generate_cdrs(
                f"Target binding site: {target_motif}",
                vh['FR1'] + "X" * 10 + vh['FR2'],  # Use template for conditioning
                3
            )
            light_cdrs = self._generate_cdrs(
                f"Light chain CDRs for {target_motif}",
                vl['FR1'] + "X" * 8 + vl['FR2'],
                3
            )
            
            if len(heavy_cdrs) >= 3 and len(light_cdrs) >= 3:
                # Assemble full sequence
                sequence = self._assemble_antibody(heavy_cdrs[:3], light_cdrs[:3])
                validation_score = self._validate_sequence(sequence)
                
                # Only include sequences that pass validation threshold
                if validation_score >= 0.7:  # Higher threshold for quality
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
