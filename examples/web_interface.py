"""
Web interface example for Healdette antibody sequence generation pipeline.
"""
from flask import Flask, render_template, request, jsonify
from modules.weighted_validator import WeightedSequenceValidator
from modules.config_validator import ConfigValidator
from modules.sequence_generator import SequenceGenerator
from flask_wtf.csrf import CSRFProtect
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_talisman import Talisman

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your-secret-key-here'  # Change this in production
csrf = CSRFProtect(app)
limiter = Limiter(app=app, key_func=get_remote_address)
Talisman(app)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/generate', methods=['POST'])
@limiter.limit("10 per minute")  # Rate limiting for security
def generate_sequences():
    try:
        # Get configuration from request
        config = request.get_json()
        
        # Validate configuration
        config_validator = ConfigValidator()
        validation_result = config_validator.validate_config(config)
        
        if not validation_result['valid']:
            return jsonify({
                'error': 'Invalid configuration',
                'details': validation_result['errors']
            }), 400
        
        # Generate sequences
        generator = SequenceGenerator(config)
        num_sequences = request.args.get('num_sequences', default=5, type=int)
        sequences = generator.generate_sequences(num_sequences=num_sequences)
        
        # Validate sequences
        validator = WeightedSequenceValidator(config)
        results = []
        
        for seq in sequences:
            validation = validator.validate_sequence(seq)
            if validation['valid']:
                results.append({
                    'sequence': seq,
                    'validation': validation
                })
        
        return jsonify({
            'sequences': results
        })
        
    except Exception as e:
        return jsonify({
            'error': 'Internal server error',
            'details': str(e)
        }), 500

if __name__ == '__main__':
    app.run(debug=False)  # Set debug=False in production