"""Flask Web Application - 11B NMR Chemical Shift Prediction System"""

import os
import sys
from flask import Flask, render_template, request, jsonify, send_file, send_from_directory
from flask_cors import CORS
from datetime import datetime
import io
import csv

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(__file__))

from config import DevelopmentConfig
from core.predictor import BoronNMRPredictor
from database.models import init_db, save_prediction, get_history, get_prediction_by_id
from utils.validators import validate_input
from utils.exceptions import InvalidSMILESError, ValidationError, PredictionError

# Create Flask application
app = Flask(__name__, template_folder='templates', static_folder='static')
app.config.from_object(DevelopmentConfig)

# Enable CORS
CORS(app)

# Global predictor instance (initialized at application startup)
predictor = None
_initialized = False


def initialize_app():
    """Application initialization"""
    global predictor, _initialized

    if _initialized:
        return

    _initialized = True

    print("\n" + "="*60)
    print("11B NMR Prediction System - Web UI")
    print("="*60)

    # Initialize database
    try:
        init_db(app.config['DATABASE_PATH'])
    except Exception as e:
        print(f"✗ Database initialization failed: {e}")
        return

    # Create image output directory
    os.makedirs(app.config['IMAGE_DIR'], exist_ok=True)

    # Load model
    print("\nLoading model...")
    try:
        predictor = BoronNMRPredictor(
            model_dir=app.config['MODEL_DIR'],
            device=app.config['DEVICE'],
            hidden_dim=app.config['HIDDEN_DIM'],
            dropout=app.config['DROPOUT'],
            solvent_dim=app.config['SOLVENT_DIM'],
            ml_feature_dim=app.config['ML_FEATURE_DIM'],
            ml_hidden_dim=app.config['ML_HIDDEN_DIM']
        )
        print("✓ Model loaded successfully\n")
    except Exception as e:
        print(f"✗ Model loading failed: {e}\n")
        raise


# Flask 3.0+ compatible initialization
@app.before_request
def before_request():
    """Execute before each request"""
    initialize_app()


# ============================================================================
# Routes - Pages
# ============================================================================

@app.route('/')
def index():
    """Main page"""
    solvents = app.config['SUPPORTED_SOLVENTS']
    return render_template('index.html', solvents=solvents)


@app.route('/history')
def history_page():
    """Prediction history page"""
    return render_template('history.html')


# ============================================================================
# Ketcher Editor Routes
# ============================================================================

@app.route('/ketcher')
def ketcher_standalone():
    """Ketcher Standalone editor - returns modified index.html"""
    ketcher_path = os.path.join(
        os.path.dirname(__file__),
        'static/lib/ketcher/ketcher-standalone'
    )

    # Read the original index.html
    index_path = os.path.join(ketcher_path, 'index.html')
    with open(index_path, 'r', encoding='utf-8') as f:
        html = f.read()

    # The Ketcher iframe is loaded at .../ketcher path
    # Relative paths ./static/ inside the iframe need to be mapped to Flask's /ketcher/static/
    # Since the iframe URL is .../ketcher (no trailing slash), the browser resolves ./xxx as .../xxx
    # Solution: inject a <base> tag into the iframe HTML <head> to set the correct base path
    base_inject = '<head><base href="ketcher/">'
    html = html.replace('<head>', base_inject, 1)

    return html


@app.route('/ketcher/<path:filename>')
def ketcher_static(filename):
    """Static files and assets for Ketcher Standalone"""
    ketcher_path = os.path.join(
        os.path.dirname(__file__),
        'static/lib/ketcher/ketcher-standalone'
    )
    return send_from_directory(ketcher_path, filename)


# ============================================================================
# API Routes - Prediction
# ============================================================================

@app.route('/api/predict', methods=['POST'])
def api_predict():
    """
    Prediction API endpoint

    Request JSON:
    {
        "molecule_smiles": "OB(O)c1ccccc1",
        "solvent": "CDCl3"
    }

    Response JSON:
    {
        "success": true,
        "prediction_id": 1,
        "canonical_smiles": "OB(O)c1ccccc1",
        "predictions": [
            {"atom_index": 1, "element": "B", "ppm": 28.53}
        ],
        "image_url": "/api/image/prediction_abc123.png",
        "num_borons": 1,
        "timestamp": "2025-12-24T10:30:00"
    }
    """
    try:
        if not predictor:
            return jsonify({
                'success': False,
                'error': 'Model not loaded, please try again later'
            }), 500

        # 1. Parse request data
        data = request.get_json()
        if not data:
            return jsonify({
                'success': False,
                'error': 'Please provide JSON data'
            }), 400

        mol_smiles = data.get('molecule_smiles', '').strip()
        solvent_name = data.get('solvent', '').strip()

        # 2. Validate inputs
        if not mol_smiles:
            return jsonify({
                'success': False,
                'error': 'molecule_smiles cannot be empty'
            }), 400

        if not solvent_name:
            return jsonify({
                'success': False,
                'error': 'solvent cannot be empty'
            }), 400

        try:
            canonical_smiles, solvent_smiles = validate_input(
                mol_smiles, solvent_name, app.config['SUPPORTED_SOLVENTS']
            )
        except (InvalidSMILESError, ValidationError) as e:
            return jsonify({
                'success': False,
                'error': str(e)
            }), 400

        # 3. Call predictor (V3 model uses solvent name, not SMILES)
        try:
            result = predictor.predict(canonical_smiles, solvent_name)
        except PredictionError as e:
            return jsonify({
                'success': False,
                'error': f'Prediction failed: {str(e)}'
            }), 500

        # 4. Generate molecule image
        try:
            image_path = predictor.generate_molecule_image(
                result['mol_object'],
                [p['atom_index'] for p in result['predictions']],
                result['predictions'],
                app.config['IMAGE_DIR']
            )
            image_filename = os.path.basename(image_path)
        except Exception as e:
            app.logger.warning(f"Molecule image generation failed: {e}")
            image_filename = None

        # 5. Save to database
        try:
            prediction_id = save_prediction(
                db_path=app.config['DATABASE_PATH'],
                mol_smiles=canonical_smiles,
                solvent_name=solvent_name,
                predictions=result['predictions'],
                image_path=image_filename if image_filename else ''
            )
        except Exception as e:
            app.logger.warning(f"Database save failed: {e}")
            prediction_id = None

        # 6. Return result
        return jsonify({
            'success': True,
            'prediction_id': prediction_id,
            'canonical_smiles': canonical_smiles,
            'predictions': result['predictions'],
            'image_url': f"/api/image/{image_filename}" if image_filename else None,
            'num_borons': result['num_borons'],
            'timestamp': datetime.now().isoformat()
        })

    except Exception as e:
        app.logger.error(f"Prediction error: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Internal server error: {str(e)}'
        }), 500


# ============================================================================
# API Routes - History
# ============================================================================

@app.route('/api/history')
def api_history():
    """
    Get prediction history API

    Query parameters:
    - limit: maximum number of records to return (default 50)

    Response JSON:
    {
        "success": true,
        "records": [
            {
                "id": 1,
                "timestamp": "2025-12-24T10:30:00",
                "mol_smiles": "OB(O)c1ccccc1",
                "solvent_name": "CDCl3",
                "predictions": [{"atom_index": 1, "element": "B", "ppm": 28.53}],
                "image_path": "prediction_abc123.png"
            }
        ]
    }
    """
    try:
        limit = request.args.get('limit', app.config['HISTORY_LIMIT'], type=int)
        limit = min(limit, 200)  # return at most 200 records

        records = get_history(app.config['DATABASE_PATH'], limit=limit)

        return jsonify({
            'success': True,
            'records': records
        })

    except Exception as e:
        app.logger.error(f"Failed to retrieve history: {e}")
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve history: {str(e)}'
        }), 500


@app.route('/api/prediction/<int:prediction_id>')
def api_get_prediction(prediction_id):
    """Get a single prediction record by ID"""
    try:
        record = get_prediction_by_id(app.config['DATABASE_PATH'], prediction_id)
        return jsonify({
            'success': True,
            'record': record
        })

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 404


# ============================================================================
# API Routes - File Serving
# ============================================================================

@app.route('/api/image/<filename>')
def api_serve_image(filename):
    """Serve molecule structure images"""
    try:
        # Security check: prevent path traversal
        if '..' in filename or '/' in filename:
            return jsonify({'error': 'Invalid filename'}), 400

        img_path = os.path.join(app.config['IMAGE_DIR'], filename)

        if not os.path.exists(img_path):
            return jsonify({'error': 'Image not found'}), 404

        return send_file(img_path, mimetype='image/png')

    except Exception as e:
        app.logger.error(f"Image serving failed: {e}")
        return jsonify({'error': 'File serving error'}), 500


# ============================================================================
# API Routes - Downloads
# ============================================================================

@app.route('/api/download/csv/<int:prediction_id>')
def api_download_csv(prediction_id):
    """Download prediction results as CSV"""
    try:
        record = get_prediction_by_id(app.config['DATABASE_PATH'], prediction_id)

        # Generate CSV
        output = io.StringIO()
        writer = csv.writer(output)

        # Write header
        writer.writerow(['Atom Index', 'Element', 'Chemical Shift (ppm)'])

        # Write data
        for pred in record['predictions']:
            writer.writerow([
                pred['atom_index'],
                pred['element'],
                f"{pred['ppm']:.2f}"
            ])

        # Write metadata
        writer.writerow([])
        writer.writerow(['Canonical SMILES', record['mol_smiles']])
        writer.writerow(['Solvent', record['solvent_name']])
        writer.writerow(['Timestamp', record['timestamp']])

        # Return file
        output.seek(0)
        return send_file(
            io.BytesIO(output.getvalue().encode('utf-8')),
            mimetype='text/csv',
            as_attachment=True,
            download_name=f"prediction_{prediction_id}.csv"
        )

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/download/json/<int:prediction_id>')
def api_download_json(prediction_id):
    """Download prediction results as JSON"""
    try:
        record = get_prediction_by_id(app.config['DATABASE_PATH'], prediction_id)

        return send_file(
            io.BytesIO(
                jsonify(record).get_data(as_text=False).encode('utf-8')
            ),
            mimetype='application/json',
            as_attachment=True,
            download_name=f"prediction_{prediction_id}.json"
        )

    except Exception as e:
        return jsonify({'error': str(e)}), 500


# ============================================================================
# Error Handlers
# ============================================================================

@app.errorhandler(404)
def not_found(error):
    """404 error handler"""
    return jsonify({'error': 'Page not found'}), 404


@app.errorhandler(500)
def internal_error(error):
    """500 error handler"""
    app.logger.error(f"Internal error: {error}")
    return jsonify({'error': 'Internal server error'}), 500


@app.errorhandler(405)
def method_not_allowed(error):
    """405 method not allowed error handler"""
    return jsonify({'error': 'Method not allowed'}), 405


# ============================================================================
# Application Entry Point
# ============================================================================

if __name__ == '__main__':
    # Ensure initialization runs within the application context
    with app.app_context():
        initialize_app()

    # Start Flask application
    app.run(
        host='0.0.0.0',
        port=5000,
        debug=app.config['DEBUG'],
        use_reloader=True
    )
