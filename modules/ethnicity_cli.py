import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Any
from ethnicity_manager import EthnicityManager

def parse_ethnicity_weights(weights_str: str) -> Dict[str, float]:
    """Parse ethnicity weights from string format 'ethnicity1:weight1,ethnicity2:weight2'"""
    weights = {}
    pairs = weights_str.split(',')
    for pair in pairs:
        try:
            ethnicity, weight = pair.split(':')
            weights[ethnicity.strip()] = float(weight.strip())
        except ValueError:
            print(f"Error: Invalid weight format in '{pair}'. Expected format: ethnicity:weight")
            sys.exit(1)
    return weights

def main():
    parser = argparse.ArgumentParser(description='Ethnicity Configuration Manager CLI')
    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Upload command
    upload_parser = subparsers.add_parser('upload', help='Upload a new ethnicity configuration')
    upload_parser.add_argument('name', help='Name of the ethnicity')
    upload_parser.add_argument('config', type=argparse.FileType('r'), help='Configuration file path')

    # Generate command
    generate_parser = subparsers.add_parser('generate', help='Generate a merged configuration')
    generate_parser.add_argument('weights', help='Ethnicity weights (format: ethnicity1:weight1,ethnicity2:weight2)')
    generate_parser.add_argument('output', help='Output file path')

    # List command
    subparsers.add_parser('list', help='List available ethnicity configurations')

    args = parser.parse_args()
    
    manager = EthnicityManager()

    if args.command == 'upload':
        try:
            config = json.load(args.config)
            if manager.validate_config(config):
                saved_path = manager.save_config(config, args.name)
                print(f"Configuration saved to {saved_path}")
            else:
                print("Error: Invalid configuration format")
                sys.exit(1)
        except json.JSONDecodeError:
            print("Error: Invalid JSON format")
            sys.exit(1)

    elif args.command == 'generate':
        try:
            weights = parse_ethnicity_weights(args.weights)
            merged_config = manager.merge_configs(weights)
            with open(args.output, 'w') as f:
                json.dump(merged_config, f, indent=2)
            print(f"Generated configuration saved to {args.output}")
        except Exception as e:
            print(f"Error generating configuration: {str(e)}")
            sys.exit(1)

    elif args.command == 'list':
        configs = list(Path(manager.config_dir).glob("*.json"))
        if configs:
            print("Available ethnicity configurations:")
            for config in configs:
                print(f"  - {config.stem}")
        else:
            print("No ethnicity configurations found")

if __name__ == '__main__':
    main()