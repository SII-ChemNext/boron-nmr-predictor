# Classical Machine Learning Benchmarks

This directory contains the configuration files and results for the classical machine learning baseline experiments. All experiments were conducted using **[Chemia](https://github.com/flyben97/Chemia)**, an open-source framework for molecular property prediction.

## Framework: Chemia

Chemia provides a YAML-driven workflow for training, evaluating, and comparing classical ML models on molecular datasets. It handles molecular featurization, data splitting, Optuna-based hyperparameter optimization, and result logging automatically.

## Installation

```bash
# 1. Clone the Chemia repository
git clone https://github.com/flyben97/Chemia.git
cd Chemia

# 2. Create and activate the conda environment
conda env create -f environment.yml
conda activate chemia

# 3. Install the package
pip install -e .
```

## Usage

All experiments are controlled through YAML configuration files. To reproduce any experiment:

```bash
# From the Chemia root directory, with the chemia environment activated
conda activate chemia
python scripts/run_training_only.py /path/to/config.yaml
```

### Example

```bash
python scripts/run_training_only.py configs/morgan_rdkit_seed42.yaml
```

## Dataset

The experiments use a **subset of the full dataset** containing only compounds with chemically equivalent boron atoms (i.e., molecules where all boron sites share the same chemical environment and thus produce a single unambiguous ¹¹B NMR shift). This subset comprises **5,217 entries** across 11 solvents, compared to 5,490 entries in the complete dataset.

The data file is located at `data/data_b11_all_single_ppm.csv` with three columns:

| Column | Description |
|---|---|
| `Smiles` | SMILES string of the boron-containing compound |
| `solvent` | SMILES string of the deuterated solvent |
| `ppm_values` | ¹¹B NMR chemical shift (ppm) |

Both `Smiles` and `solvent` are independently featurized and concatenated as model input.

## Configuration Files

The `configs/` directory contains 33 configuration files covering **11 molecular representations × 3 random seeds** (seed 42, 43, 44). Each seed run uses an independent 80/10/10 train/validation/test split.

| Config file prefix | Molecular representation | Feature type |
|---|---|---|
| `maccs_*` | MACCS Keys | 167-bit structural keys |
| `rdkit_descriptors_*` | RDKit Descriptors | 200 physicochemical descriptors |
| `molt5_*` | MolT5 Embeddings | Text-based molecular embeddings |
| `unimol_*` | UniMol v2 (310M) | 3D-geometry-aware embeddings |
| `morgan_*` | Morgan Fingerprints | Circular fingerprints (radius=2, 1024 bits) |
| `maccs_rdkit_*` | MACCS + RDKit | Concatenated fingerprints + descriptors |
| `molt5_rdkit_*` | MolT5 + RDKit | Text embeddings + descriptors |
| `molt5_morgan_*` | MolT5 + Morgan | Text embeddings + fingerprints |
| `molt5_maccs_*` | MolT5 + MACCS | Text embeddings + structural keys |
| `unimol_rdkit_*` | UniMol + RDKit | 3D embeddings + descriptors |
| `morgan_rdkit_*` | Morgan + RDKit | Fingerprints + descriptors **(best)** |

Both the solute SMILES and solvent SMILES are featurized independently using the same representation and the resulting features are concatenated.

### Algorithms

Each configuration trains and evaluates **7 regression algorithms** simultaneously:

- XGBoost
- Random Forest
- Ridge Regression
- CatBoost
- HistGradientBoosting
- LightGBM
- Extra Trees

Hyperparameters are optimized using [Optuna](https://optuna.org/) with the TPE sampler (30 trials, MedianPruner).

### Key Configuration Parameters

```yaml
training:
  n_trials: 30
  optuna_config:
    sampler: "TPE"
    pruner: "MedianPruner"
    n_startup_trials: 10

split_mode: "train_valid_test"
split_config:
  train_valid_test:
    valid_size: 0.1
    test_size: 0.1
    random_state: 42   # changed to 43/44 for other seeds

evaluation:
  primary_metric: "r2"
  additional_metrics: ["rmse", "mae", "mape"]
```

## Results

Full results for all 11 representations × 7 algorithms × 3 seeds are in `results/model_comparison.csv`.

The best-performing combination is **Morgan + RDKit Descriptors with XGBoost** (seed 44), whose trained model weights and data splits are saved in `Morgan-RDKit_seed44_xgboost/`.

| Representation | Best Algorithm | Mean Test R² | Mean Test RMSE (ppm) | Mean Test MAE (ppm) |
|---|---|---|---|---|
| Morgan + RDKit | XGBoost | **0.9212** | **6.18** | **2.76** |
| MACCS | XGBoost | 0.9107 | 6.60 | 2.91 |
| UniMol + RDKit | HistGradientBoosting | 0.8985 | 7.02 | 3.39 |
| RDKit Descriptors | XGBoost | 0.8974 | 7.06 | 3.16 |
| MACCS + RDKit | XGBoost | 0.8933 | 6.91 | 3.64 |
| Morgan | XGBoost | 0.8820 | 7.61 | 3.63 |
| MolT5 + RDKit | HistGradientBoosting | 0.8772 | 7.73 | 3.89 |
| MolT5 + Morgan | LightGBM | 0.8751 | 7.78 | 4.06 |
| MolT5 + MACCS | XGBoost | 0.8712 | 7.91 | 3.83 |
| MolT5 | CatBoost | 0.7601 | 10.84 | 6.46 |
| UniMol | LightGBM | 0.3698 | 16.54 | 11.73 |

*Each row shows the best algorithm selected by mean Test R² across seeds 42, 43, and 44.*
*The trained model saved in `Morgan-RDKit_seed44_xgboost/` corresponds to the best-performing seed of the top configuration.*

## SHAP Interpretability Analysis

SHAP analysis was performed on the best-performing model (Morgan + RDKit / XGBoost, seed 44) using `shap.TreeExplainer` to identify the most influential features for ¹¹B NMR chemical shift prediction. All code is in `Morgan-RDKit_seed44_xgboost/SHAP/`.

> **Note:** The `data_splits/` directory is not included in this repository due to file size. To regenerate it, run Chemia with `configs/morgan_rdkit_seed44.yaml` first. The featurized X matrices will be saved automatically to `Morgan-RDKit_seed44_xgboost/data_splits/`.

```bash
cd Morgan-RDKit_seed44_xgboost/SHAP/
python run_shap_analysis.py            # global SHAP values + beeswarm summary plot
python run_morgan_analysis.py          # Morgan bit substructure analysis
python analyze_mrlow_hybridization.py  # BCUT2D_MRLOW vs. boron hybridization
```

## Reproducing All Experiments

To reproduce the full benchmark (11 representations × 3 seeds = 33 runs):

```bash
conda activate chemia
cd Chemia/

for config in /path/to/classical_ml/configs/*.yaml; do
    python scripts/run_training_only.py "$config"
done
```

## Citation

If you use Chemia in your work, please cite:

```
@software{chemia,
  author = {flyben97},
  title  = {Chemia: A Framework for Molecular Property Prediction},
  url    = {https://github.com/flyben97/Chemia},
}
```
