import json
import os
from pathlib import Path
from typing import Dict, Any, Optional

class EthnicityManager:
    def __init__(self, config_dir: str = "configs"):
        self.config_dir = Path(config_dir)
        self.config_dir.mkdir(exist_ok=True)
        
    def load_config(self, config_path: str) -> Dict[str, Any]:
        """Load an ethnicity configuration from a file."""
        with open(config_path, 'r') as f:
            return json.load(f)
    
    def save_config(self, config: Dict[str, Any], ethnicity_name: str) -> str:
        """Save an ethnicity configuration to a file."""
        config_path = self.config_dir / f"{ethnicity_name.lower().replace(' ', '_')}.json"
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)
        return str(config_path)
    
    def validate_config(self, config: Dict[str, Any]) -> bool:
        """Validate that a configuration has all required fields."""
        required_fields = {
            "population_weights",
            "binding_parameters",
            "validation_thresholds"
        }
        return all(field in config for field in required_fields)
    
    def merge_configs(self, configs: Dict[str, float]) -> Dict[str, Any]:
        """Merge multiple ethnicity configs with given weights."""
        merged = {
            "population_weights": {},
            "binding_parameters": {},
            "validation_thresholds": {}
        }
        
        total_weight = sum(configs.values())
        normalized_weights = {k: v/total_weight for k, v in configs.items()}
        
        for config_path, weight in configs.items():
            config = self.load_config(config_path)
            if not self.validate_config(config):
                raise ValueError(f"Invalid configuration in {config_path}")
                
            # Merge each section with weights
            for section in merged.keys():
                for key, value in config[section].items():
                    if key not in merged[section]:
                        merged[section][key] = 0
                    merged[section][key] += value * normalized_weights[config_path]
        
        return merged