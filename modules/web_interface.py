from flask import Flask, request, render_template, jsonify
from flask_wtf.csrf import CSRFProtect
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_talisman import Talisman
import os
import json
from pathlib import Path
from ethnicity_manager import EthnicityManager
from .security import SecurityValidator
from .logging_manager import SecurityLogger
from .sanitizer import OutputSanitizer
from .backup_manager import BackupManager

# Initialize components
security_logger = SecurityLogger()
sanitizer = OutputSanitizer()
backup_manager = BackupManager('configs')

def get_client_info():
    """Get client information for logging."""
    return {
        'ip': request.remote_addr,
        'user_agent': request.user_agent.string,
        'method': request.method,
        'path': request.path
    }

# Security configuration
CSP = {
    'default-src': "'self'",
    'script-src': "'self' 'unsafe-inline'",
    'style-src': "'self' 'unsafe-inline'",
}

app = Flask(__name__)
app.config['SECRET_KEY'] = os.urandom(32)
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024  # 1MB max file size

# Security middleware
csrf = CSRFProtect(app)
Talisman(app, content_security_policy=CSP, force_https=True)
limiter = Limiter(
    app=app,
    key_func=get_remote_address,
    default_limits=["100 per day", "10 per minute"]
)

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
    client_info = get_client_info()
    try:
        data = request.get_json()
        if not data or 'weights' not in data:
            security_logger.log_validation_failure(
                'CONFIG_GENERATION',
                data,
                'Missing weights in request'
            )
            return jsonify({'error': 'No weights provided in request'}), 400
        
        weights = data['weights']
        if not weights or not isinstance(weights, dict):
            security_logger.log_validation_failure(
                'CONFIG_GENERATION',
                weights,
                'Invalid weights format'
            )
            return jsonify({'error': 'Invalid weights format'}), 400
        
        # Convert paths to absolute paths in the configs directory
        config_weights = {
            str(Path(manager.config_dir) / f"{ethnicity}.json"): weight
            for ethnicity, weight in weights.items()
        }
        
        security_logger.log_audit_event(
            action='GENERATE_CONFIG',
            user=request.remote_addr,
            resource='config_generation',
            details={'weights': weights}
        )
        
        merged_config = manager.merge_configs(config_weights)
        output_path = Path('generated_config.json')
        with open(output_path, 'w') as f:
            json.dump(merged_config, f, indent=2)
            
        security_logger.log_access_attempt(
            resource='config_generation',
            access_type='WRITE',
            status='SUCCESS',
            client_info=client_info
        )
        
        return jsonify({
            'message': f'Configuration generated at {output_path}',
            'config': merged_config
        })
    except json.JSONDecodeError:
        security_logger.log_security_event(
            'INVALID_INPUT',
            'Invalid JSON in request',
            client_info,
            'ERROR'
        )
        return jsonify({'error': 'Invalid JSON in request'}), 400
    except ValueError as e:
        security_logger.log_security_event(
            'VALIDATION_ERROR',
            str(e),
            client_info,
            'ERROR'
        )
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        security_logger.log_security_event(
            'SYSTEM_ERROR',
            str(e),
            client_info,
            'ERROR'
        )
        return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    app.run(debug=True)