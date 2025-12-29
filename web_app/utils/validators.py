"""输入验证工具"""

from rdkit import Chem
from utils.exceptions import InvalidSMILESError, ValidationError


def validate_smiles(smiles):
    """
    验证 SMILES 合法性

    Args:
        smiles (str): SMILES 字符串

    Returns:
        bool: 验证通过返回 True

    Raises:
        InvalidSMILESError: SMILES 无效
    """
    if not smiles or not isinstance(smiles, str):
        raise InvalidSMILESError("SMILES 不能为空")

    smiles = smiles.strip()
    if not smiles:
        raise InvalidSMILESError("SMILES 不能为空")

    # 尝试解析 SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"无法解析的 SMILES: {smiles}")

    # 检查是否包含硼原子
    has_boron = any(atom.GetSymbol() == 'B' for atom in mol.GetAtoms())
    if not has_boron:
        raise InvalidSMILESError("分子中必须包含至少一个硼原子 (B)")

    return True


def validate_solvent(solvent_name, supported_solvents):
    """
    验证溶剂并返回对应的 SMILES

    Args:
        solvent_name (str): 溶剂名称
        supported_solvents (dict): 支持的溶剂字典

    Returns:
        str: 溶剂的 SMILES 字符串

    Raises:
        ValidationError: 溶剂不支持
    """
    if solvent_name not in supported_solvents:
        raise ValidationError(
            f"不支持的溶剂: {solvent_name}. "
            f"支持的溶剂有: {', '.join(supported_solvents.keys())}"
        )

    return supported_solvents[solvent_name]


def validate_input(smiles, solvent_name, supported_solvents):
    """
    验证预测输入

    Args:
        smiles (str): 分子 SMILES
        solvent_name (str): 溶剂名称
        supported_solvents (dict): 支持的溶剂字典

    Returns:
        tuple: (canonical_smiles, solvent_smiles)

    Raises:
        InvalidSMILESError: SMILES 无效
        ValidationError: 溶剂无效
    """
    # 验证 SMILES
    validate_smiles(smiles)

    # 标准化 SMILES
    mol = Chem.MolFromSmiles(smiles)
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    # 验证溶剂
    solvent_smiles = validate_solvent(solvent_name, supported_solvents)

    return canonical_smiles, solvent_smiles
