"""
Celtic/British ancestry configuration for Healdette pipeline.
"""

# HLA Profiles for different Celtic/British ancestries
HLA_PROFILES = {
    "celtic_british": {
        "hla_a": ["A*01:01", "A*03:01", "A*02:01"],  # Most common Celtic alleles
        "hla_b": ["B*07:02", "B*08:01", "B*44:02"],  # British Isles specific
        "hla_c": ["C*07:01", "C*04:01", "C*06:02"]   # Celtic markers
    },
    "irish_dominant": {
        "hla_a": ["A*01:01", "A*02:01", "A*03:01"],  # Irish frequencies
        "hla_b": ["B*08:01", "B*07:02", "B*44:02"],  # Irish-specific
        "hla_c": ["C*07:01", "C*04:01", "C*06:02"]   # Irish markers
    },
    "scottish": {
        "hla_a": ["A*02:01", "A*01:01", "A*03:01"],  # Scottish frequencies
        "hla_b": ["B*07:02", "B*44:02", "B*08:01"],  # Highland markers
        "hla_c": ["C*07:01", "C*04:01", "C*05:01"]   # Scottish common
    }
}

# Celtic-specific binding motifs and weights
CELTIC_MOTIFS = {
    "RF": 0.15,  # Strong Celtic association
    "KW": 0.12,  # British Isles marker
    "WY": 0.10,  # Common in Celtic populations
    "YF": 0.08   # Highland Scottish marker
}

# Population-specific biokinetic factors
BIOKINETIC_FACTORS = {
    "celtic_british": {
        "metabolism": 1.0,      # Baseline
        "clearance": 1.0,       # Standard
        "absorption": 1.0       # Reference
    },
    "irish_dominant": {
        "metabolism": 1.02,     # Slightly enhanced
        "clearance": 1.05,      # Faster clearance
        "absorption": 1.03      # Better absorption
    },
    "scottish": {
        "metabolism": 0.98,     # Slightly reduced
        "clearance": 0.95,      # Slower clearance
        "absorption": 0.97      # Standard absorption
    }
}

# Ancestry-specific immune response modifiers
IMMUNE_MODIFIERS = {
    "celtic_british": {
        "baseline": 1.0,
        "inflammation": 1.02,
        "antibody_production": 1.05
    },
    "irish_dominant": {
        "baseline": 1.05,
        "inflammation": 1.08,
        "antibody_production": 1.10
    },
    "scottish": {
        "baseline": 1.02,
        "inflammation": 1.05,
        "antibody_production": 1.07
    }
}

# Validation thresholds for Celtic populations
VALIDATION_THRESHOLDS = {
    "min_celtic_score": 0.2,
    "min_population_coverage": 0.6,
    "min_hla_matches": 3,
    "required_motifs": 1
}

# Population-specific disease susceptibilities
DISEASE_SUSCEPTIBILITIES = {
    "celtic_british": {
        "autoimmune": "high",
        "cardiovascular": "moderate",
        "metabolic": "standard"
    },
    "irish_dominant": {
        "autoimmune": "very_high",
        "cardiovascular": "high",
        "metabolic": "moderate"
    },
    "scottish": {
        "autoimmune": "high",
        "cardiovascular": "high",
        "metabolic": "moderate"
    }
}

def get_celtic_config(ancestry_profile=None):
    """
    Get configuration settings for specified Celtic ancestry profile.
    
    Args:
        ancestry_profile: List of ancestry tags (e.g., ["celtic_british", "irish_dominant"])
        
    Returns:
        Dict containing configuration settings
    """
    if ancestry_profile is None:
        ancestry_profile = ["celtic_british"]
        
    primary_ancestry = ancestry_profile[0]
    
    config = {
        "hla_profile": HLA_PROFILES.get(primary_ancestry, HLA_PROFILES["celtic_british"]),
        "binding_motifs": CELTIC_MOTIFS,
        "biokinetics": BIOKINETIC_FACTORS.get(primary_ancestry, BIOKINETIC_FACTORS["celtic_british"]),
        "immune_modifiers": IMMUNE_MODIFIERS.get(primary_ancestry, IMMUNE_MODIFIERS["celtic_british"]),
        "validation": VALIDATION_THRESHOLDS,
        "disease_risks": DISEASE_SUSCEPTIBILITIES.get(primary_ancestry, 
                                                    DISEASE_SUSCEPTIBILITIES["celtic_british"])
    }
    
    # Adjust for mixed ancestry if multiple profiles
    if len(ancestry_profile) > 1:
        secondary_ancestry = ancestry_profile[1]
        primary_weight = 0.7
        secondary_weight = 0.3
        
        # Blend biokinetic factors
        for key in config["biokinetics"]:
            primary_value = BIOKINETIC_FACTORS[primary_ancestry][key]
            secondary_value = BIOKINETIC_FACTORS[secondary_ancestry][key]
            config["biokinetics"][key] = (primary_value * primary_weight + 
                                        secondary_value * secondary_weight)
        
        # Blend immune modifiers
        for key in config["immune_modifiers"]:
            primary_value = IMMUNE_MODIFIERS[primary_ancestry][key]
            secondary_value = IMMUNE_MODIFIERS[secondary_ancestry][key]
            config["immune_modifiers"][key] = (primary_value * primary_weight + 
                                             secondary_value * secondary_weight)
    
    return config