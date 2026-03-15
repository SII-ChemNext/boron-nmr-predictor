#!/bin/bash

# 11B NMR Prediction System - Web UI startup script

set -e

echo "=================================================="
echo "11B NMR Prediction System - Web UI"
echo "=================================================="
echo ""

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "Working directory: $SCRIPT_DIR"
cd "$SCRIPT_DIR"

# Check Python version
echo "Checking Python version..."
if ! command -v python3 &> /dev/null; then
    echo "✗ Error: Python 3 not found"
    exit 1
fi

PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "✓ Python version: $PYTHON_VERSION"
echo ""

# Check virtual environment
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
    echo "✓ Virtual environment created"
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate
echo "✓ Virtual environment activated"
echo ""

# Check dependencies
echo "Checking dependencies..."
if [ ! -f "requirements.txt" ]; then
    echo "✗ Error: requirements.txt not found"
    exit 1
fi

pip install -q -r requirements.txt
echo "✓ Dependencies installed"
echo ""

# Check model files
echo "Checking model files..."
MISSING_MODELS=0
for i in {1..5}; do
    if [ ! -f "models/model_v3_fold_$i.pth" ]; then
        echo "✗ Missing model file: models/model_v3_fold_$i.pth"
        MISSING_MODELS=1
    fi
done

if [ -f "core/model.py" ] && [ -f "core/features.py" ]; then
    echo "✓ Model file check complete"
else
    echo "✗ Error: required model source files are missing"
    exit 1
fi
echo ""

# Create required directories
echo "Creating required directories..."
mkdir -p database static/img
echo "✓ Directories created"
echo ""

# Check database
if [ ! -f "database/predictions.db" ]; then
    echo "Initializing database..."
    python3 << EOF
from database.models import init_db
init_db('database/predictions.db')
EOF
    echo "✓ Database initialized"
else
    echo "✓ Database already exists"
fi
echo ""

# Start application
echo "=================================================="
echo "Starting application..."
echo "=================================================="
echo ""
echo "Service URL: http://localhost:5000"
echo ""
echo "Press Ctrl+C to stop"
echo ""

# Check whether to use Gunicorn (production mode)
if [ "$1" == "production" ]; then
    echo "Starting in production mode (Gunicorn)..."
    pip install -q gunicorn
    gunicorn -w 4 -b 0.0.0.0:5000 --timeout 120 app:app
else
    echo "Starting in development mode (Flask)..."
    python3 app.py
fi
