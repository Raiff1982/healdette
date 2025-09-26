from flask import Flask, request, render_template, jsonify
import os
import json
from pathlib import Path
from ethnicity_manager import EthnicityManager

app = Flask(__name__)
manager = EthnicityManager()

@app.route('/')
def index():
    configs = list(Path(manager.config_dir).glob("*.json"))
    ethnicities = [config.stem for config in configs]
    return render_template('index.html', ethnicities=ethnicities)

@app.route('/upload', methods=['POST'])
def upload_config():
    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400
    
    file = request.files['file']
    ethnicity_name = request.form.get('name')
    
    if not ethnicity_name:
        return jsonify({'error': 'No ethnicity name provided'}), 400
    
    if file:
        try:
            config = json.load(file)
            if manager.validate_config(config):
                saved_path = manager.save_config(config, ethnicity_name)
                return jsonify({'message': f'Configuration saved to {saved_path}'})
            else:
                return jsonify({'error': 'Invalid configuration format'}), 400
        except json.JSONDecodeError:
            return jsonify({'error': 'Invalid JSON format'}), 400
    
    return jsonify({'error': 'No file selected'}), 400

@app.route('/generate', methods=['POST'])
def generate_config():
    try:
        data = request.get_json()
        if not data or 'weights' not in data:
            return jsonify({'error': 'No weights provided in request'}), 400
        
        weights = data['weights']
        if not weights or not isinstance(weights, dict):
            return jsonify({'error': 'Invalid weights format'}), 400
        
        # Convert paths to absolute paths in the configs directory
        config_weights = {
            str(Path(manager.config_dir) / f"{ethnicity}.json"): weight
            for ethnicity, weight in weights.items()
        }
        
        merged_config = manager.merge_configs(config_weights)
        output_path = Path('generated_config.json')
        with open(output_path, 'w') as f:
            json.dump(merged_config, f, indent=2)
        return jsonify({
            'message': f'Configuration generated at {output_path}',
            'config': merged_config
        })
    except json.JSONDecodeError:
        return jsonify({'error': 'Invalid JSON in request'}), 400
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        app.logger.error(f'Unexpected error: {str(e)}')
        return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    app.run(debug=True)