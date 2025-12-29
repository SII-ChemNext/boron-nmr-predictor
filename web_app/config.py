import os
import torch

class Config:
    """应用配置"""

    # 基础路径
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    # 模型配置
    MODEL_DIR = os.path.join(BASE_DIR, 'models')
    DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

    # 数据库配置
    DATABASE_DIR = os.path.join(BASE_DIR, 'database')
    DATABASE_PATH = os.path.join(DATABASE_DIR, 'predictions.db')

    # 静态文件路径
    IMAGE_DIR = os.path.join(BASE_DIR, 'static', 'img')

    # Flask 配置
    DEBUG = True
    SECRET_KEY = 'boron-nmr-prediction-web-app-secret-key-2025'
    JSON_SORT_KEYS = False

    # 支持的溶剂列表（名称 -> 氘代 SMILES）
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

    # 模型参数（必须与训练时一致）
    HIDDEN_DIM = 256
    DROPOUT = 0.0472
    SOLVENT_FP_SIZE = 1024

    # 文件上传限制
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB

    # 历史记录查询限制
    HISTORY_LIMIT = 50


class DevelopmentConfig(Config):
    """开发环境配置"""
    DEBUG = True
    TESTING = False


class ProductionConfig(Config):
    """生产环境配置"""
    DEBUG = False
    TESTING = False


class TestingConfig(Config):
    """测试环境配置"""
    TESTING = True
    DEBUG = True
    DATABASE_PATH = ':memory:'  # 使用内存数据库进行测试


# 环境变量选择
def get_config():
    """根据环境变量选择配置"""
    env = os.getenv('FLASK_ENV', 'development')

    if env == 'production':
        return ProductionConfig
    elif env == 'testing':
        return TestingConfig
    else:
        return DevelopmentConfig
