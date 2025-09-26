"""
Configuration validator for Healdette's multi-ethnic antibody generation pipeline.
"""

import json
import jsonschema
from pathlib import Path
from typing import Dict, Union, Optional

class ConfigValidator:
    """Validates Healdette configuration files against the schema."""
    
    def __init__(self, schema_path: Optional[str] = None):
        """
        Initialize the validator with a schema.
        
        Args:
            schema_path: Path to the JSON schema file. If None, uses the default schema.
        """
        if schema_path is None:
            schema_path = str(Path(__file__).parent.parent / 'schemas' / 'config_schema.json')
            
        with open(schema_path, 'r') as f:
            self.schema = json.load(f)
            
    def validate(self, config: Union[str, Dict, Path]) -> Dict:
        """
        Validate a configuration against the schema.
        
        Args:
            config: Either a path to a JSON file, a dictionary, or a Path object
            
        Returns:
            Dict containing validation results
            
        Raises:
            jsonschema.exceptions.ValidationError: If validation fails
            jsonschema.exceptions.SchemaError: If the schema itself is invalid
        """
        # Load configuration if it's a file path
        if isinstance(config, (str, Path)):
            with open(config, 'r') as f:
                config_data = json.load(f)
        else:
            config_data = config
            
        # Perform validation
        try:
            jsonschema.validate(instance=config_data, schema=self.schema)
            return {
                'valid': True,
                'errors': None
            }
        except jsonschema.exceptions.ValidationError as e:
            return {
                'valid': False,
                'errors': [
                    {
                        'message': e.message,
                        'path': ' -> '.join(str(p) for p in e.path),
                        'schema_path': ' -> '.join(str(p) for p in e.schema_path)
                    }
                ]
            }
            
    def validate_population(self, population_config: Dict) -> Dict:
        """
        Validate a single population's configuration.
        
        Args:
            population_config: Dictionary containing population-specific parameters
            
        Returns:
            Dict containing validation results
        """
        population_schema = self.schema['properties']['populations']['patternProperties']['^[a-z_]+$']
        
        try:
            jsonschema.validate(instance=population_config, schema=population_schema)
            return {
                'valid': True,
                'errors': None
            }
        except jsonschema.exceptions.ValidationError as e:
            return {
                'valid': False,
                'errors': [
                    {
                        'message': e.message,
                        'path': ' -> '.join(str(p) for p in e.path),
                        'schema_path': ' -> '.join(str(p) for p in e.schema_path)
                    }
                ]
            }
            
    def get_population_template(self) -> Dict:
        """
        Get a template for a population configuration.
        
        Returns:
            Dict containing the structure for a population configuration
        """
        return {
            "ancestry_weight": 0.0,
            "binding_motifs": ["XX"],
            "biophysical_params": {
                "aromatic_content": {
                    "min": 0,
                    "max": 100
                },
                "hydrophobic_content": {
                    "min": 0,
                    "max": 100
                },
                "net_charge": {
                    "min": -20,
                    "max": 20
                }
            },
            "hla_frequencies": {
                "hla_a": {},
                "hla_b": {},
                "hla_c": {}
            },
            "notes": "Population-specific notes"
        }
        
    def validate_file(self, file_path: Union[str, Path]) -> Dict:
        """
        Validate a configuration file and provide detailed error messages.
        
        Args:
            file_path: Path to the configuration file
            
        Returns:
            Dict containing validation results with detailed messages
        """
        try:
            with open(file_path, 'r') as f:
                config_data = json.load(f)
        except json.JSONDecodeError as e:
            return {
                'valid': False,
                'errors': [
                    {
                        'message': f'Invalid JSON: {str(e)}',
                        'line': e.lineno,
                        'column': e.colno
                    }
                ]
            }
            
        return self.validate(config_data)