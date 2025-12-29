# Boron NMR Chemical Shift Predictor

English | [简体中文](README_CN.md)

<sup>11</sup>B-NMR Chemical Shift Prediction Web Application based on Graph Neural Networks

## Interface Preview

![Web Interface](docs/images/web-interface.png)

*Main interface showing the Ketcher molecular editor (left) and prediction results panel (right)*

> **Note**: If you don't see the image above, please add your screenshot to `docs/images/web-interface.png`

## About

This is a web application for predicting <sup>11</sup>B-NMR chemical shifts using deep learning. The system employs a Graph Neural Network (GNN) architecture combined with a 5-fold cross-validation ensemble model to accurately predict chemical shifts of boron-containing molecules in various deuterated solvents.

### Key Features

**Molecular Editor (Left Panel)**
- 🎨 Ketcher visual molecular editor
- ⚛️ Support for complex molecular structures
- 🔧 Complete chemical drawing tools

**Prediction Features (Right Panel)**
- 🔬 10 deuterated solvents available
- 🤖 Ensemble prediction with 5 models
- 📊 Real-time results display
- 💾 Automatic history saving

**Results Display**
- 📋 Prediction table (atom index, element, ppm value)
- 🖼️ Molecular structure image (boron atoms highlighted)
- 📥 Multi-format download (CSV, JSON, PNG)

## Tech Stack

- **Backend**: Flask 3.0
- **Deep Learning**: PyTorch 2.1 + PyTorch Geometric 2.4
- **Cheminformatics**: RDKit 2023.9
- **Frontend**: Vanilla JavaScript + Ketcher Molecular Editor
- **Database**: SQLite3

## Quick Start

### Requirements

- Python 3.8+
- 4GB+ RAM (8GB recommended)
- Modern browser (Chrome, Firefox, Safari, Edge)

### Installation

1. **Clone the repository**
```bash
git clone https://github.com/your-username/boron-nmr-predictor.git
cd boron-nmr-predictor/web_app
```

2. **Create virtual environment**
```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate  # Windows
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```

4. **Download model files**

Model files (~500MB) are not included in the repository. Download from:

- [Google Drive](https://drive.google.com/drive/folders/1HLirbH9JOf6HvgwUTIUETkpXG_U4fPYl?usp=sharing)

Place the 5 model files in `web_app/models/` directory:
```
models/
├── model_fold_1.pth
├── model_fold_2.pth
├── model_fold_3.pth
├── model_fold_4.pth
└── model_fold_5.pth
```

5. **Start the application**
```bash
python app.py
```

6. **Access the application**

Open in browser: http://localhost:5000

## Usage

### Predicting Chemical Shifts

1. Draw boron-containing molecule in the editor (or input SMILES directly)
2. Select solvent type
3. Click "Predict" button
4. View prediction results and molecular structure
5. Optional: Download results (CSV/JSON/PNG)

### Supported Solvents

| Solvent | Formula |
|---------|---------|
| CDCl3 | Deuterated Chloroform |
| C6D6 | Deuterated Benzene |
| DMSO-d6 | Deuterated DMSO |
| Acetone-d6 | Deuterated Acetone |
| CD3CN | Deuterated Acetonitrile |
| CD3OD | Deuterated Methanol |
| CD2Cl2 | Deuterated Dichloromethane |
| THF-d8 | Deuterated THF |
| Toluene-d8 | Deuterated Toluene |
| D2O | Heavy Water |

## Project Structure

```
boron-nmr-predictor/
├── README.md
├── LICENSE
├── .gitignore
└── web_app/
    ├── app.py                    # Flask main application
    ├── config.py                 # Configuration file
    ├── requirements.txt          # Python dependencies
    ├── start.sh                  # Linux/Mac startup script
    ├── start.bat                 # Windows startup script
    ├── core/                     # Core prediction logic
    │   ├── predictor.py         # Predictor class
    │   ├── model.py             # GNN model architecture
    │   └── features.py          # Feature extraction
    ├── database/                 # Database operations
    │   └── models.py
    ├── utils/                    # Utility functions
    │   ├── validators.py
    │   └── exceptions.py
    ├── templates/                # HTML templates
    │   ├── base.html
    │   ├── index.html
    │   └── history.html
    ├── static/                   # Static resources
    │   ├── css/main.css
    │   ├── js/
    │   │   ├── main.js
    │   │   ├── ketcher-integration.js
    │   │   └── result-handler.js
    │   └── img/                 # Generated molecular images
    └── models/                   # Model files (download required)
        └── model_fold_*.pth
```

## FAQ

**Q: Model loading failed?**  
A: Ensure all 5 model files are downloaded and placed in the correct location.

**Q: Ketcher editor not loading?**  
A: The app will automatically fall back to SMILES input box. Functionality is not affected.

**Q: Prediction is slow?**  
A: First prediction requires model loading (~3-5 seconds). Subsequent predictions are faster. GPU significantly accelerates prediction.

## Contributing

Issues and Pull Requests are welcome!

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this project in your research, please cite:

```bibtex
@software{boron_nmr_predictor,
  title = {Boron NMR Chemical Shift Predictor},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/your-username/boron-nmr-predictor}
}
```

## Contact

- Project Homepage: https://github.com/your-username/boron-nmr-predictor
- Issue Tracker: https://github.com/your-username/boron-nmr-predictor/issues
