import json
import os
from pathlib import Path
from typing import Dict, Any, Optional

from .security import SecurityValidator

class EthnicityManager:
    def __init__(self, config_dir: str = "configs"):
        self.config_dir = Path(config_dir).resolve()
        self.config_dir.mkdir(exist_ok=True)
        self.security = SecurityValidator()
        
    def load_config(self, config_path: str) -> Dict[str, Any]:
        """Load an ethnicity configuration from a file securely."""
        config_path = Path(config_path).resolve()
        
        if not str(config_path).startswith(str(self.config_dir)):
            raise ValueError("Access denied: File outside config directory")
            
        if not self.security.is_safe_file_size(str(config_path)):
            raise ValueError("File too large")
            
        if not self.security.get_secure_file_type(str(config_path)):
            raise ValueError("Invalid file type")
            
        with open(config_path, 'r') as f:
            config = json.load(f)
            if not self.security.is_safe_json_content(config):
                raise ValueError("Invalid JSON content")
            return config
    
    def save_config(self, config: Dict[str, Any], ethnicity_name: str) -> str:
        """Save an ethnicity configuration to a file securely."""
        if not self.security.is_safe_json_content(config):
            raise ValueError("Invalid configuration content")
            
        safe_name = self.security.sanitize_ethnicity_name(ethnicity_name)
        if not safe_name:
            raise ValueError("Invalid ethnicity name")
            
        config_path = self.security.secure_path_join(
            str(self.config_dir),
            f"{safe_name}.json"
        )
        
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)
        return config_path
    
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