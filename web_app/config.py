import os
import torch

class Config:
    """Application configuration"""

    # Base path
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    # Model configuration
    MODEL_DIR = os.path.join(BASE_DIR, 'models')
    DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

    # Database configuration
    DATABASE_DIR = os.path.join(BASE_DIR, 'database')
    DATABASE_PATH = os.path.join(DATABASE_DIR, 'predictions.db')

    # Static file path
    IMAGE_DIR = os.path.join(BASE_DIR, 'static', 'img')

    # Flask configuration
    DEBUG = True
    SECRET_KEY = 'boron-nmr-prediction-web-app-secret-key-2025'
    JSON_SORT_KEYS = False

    # Supported solvents (name -> deuterated SMILES)
    SUPPORTED_SOLVENTS = {
        'CDCl3': '[2H]C(Cl)(Cl)Cl',
        'C6D6': '[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]',
        'd6-DMSO': '[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]',
        'CD3COCD3': '[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]',
        'CD3CN': '[2H]C([2H])([2H])C#N',
        'CD3OD': '[2H]OC([2H])([2H])[2H]',
        'CD2Cl2': '[2H]C([2H])(Cl)Cl',
        'd8-THF': '[2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]',
        'd8-Toluene': '[2H]c1c([2H])c([2H])c(C([2H])([2H])[2H])c([2H])c1[2H]',
        'D2O': '[2H]O[2H]'
    }

    # Model hyperparameters (must match training configuration - V3 model)
    HIDDEN_DIM = 256
    DROPOUT = 0.012558398103042557
    SOLVENT_DIM = 32
    ML_FEATURE_DIM = 20
    ML_HIDDEN_DIM = 64

    # File upload limit
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB

    # History query limit
    HISTORY_LIMIT = 50


class DevelopmentConfig(Config):
    """Development environment configuration"""
    DEBUG = True
    TESTING = False


class ProductionConfig(Config):
    """Production environment configuration"""
    DEBUG = False
    TESTING = False


class TestingConfig(Config):
    """Testing environment configuration"""
    TESTING = True
    DEBUG = True
    DATABASE_PATH = ':memory:'  # Use in-memory database for testing


# Select configuration based on environment variable
def get_config():
    """Select configuration based on environment variable"""
    env = os.getenv('FLASK_ENV', 'development')

    if env == 'production':
        return ProductionConfig
    elif env == 'testing':
        return TestingConfig
    else:
        return DevelopmentConfig
