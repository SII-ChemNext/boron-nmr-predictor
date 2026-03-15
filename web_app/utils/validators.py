"""Input validation utilities"""

from rdkit import Chem
from utils.exceptions import InvalidSMILESError, ValidationError


def validate_smiles(smiles):
    """
    Validate SMILES string.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if validation passes

    Raises:
        InvalidSMILESError: if SMILES is invalid
    """
    if not smiles or not isinstance(smiles, str):
        raise InvalidSMILESError("SMILES cannot be empty")

    smiles = smiles.strip()
    if not smiles:
        raise InvalidSMILESError("SMILES cannot be empty")

    # Attempt to parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"Cannot parse SMILES: {smiles}")

    # Check for at least one boron atom
    has_boron = any(atom.GetSymbol() == 'B' for atom in mol.GetAtoms())
    if not has_boron:
        raise InvalidSMILESError("Molecule must contain at least one boron atom (B)")

    return True


def validate_solvent(solvent_name, supported_solvents):
    """
    Validate solvent and return its corresponding SMILES.

    Args:
        solvent_name (str): solvent name
        supported_solvents (dict): dictionary of supported solvents

    Returns:
        str: SMILES string of the solvent

    Raises:
        ValidationError: if solvent is not supported
    """
    if solvent_name not in supported_solvents:
        raise ValidationError(
            f"Unsupported solvent: {solvent_name}. "
            f"Supported solvents: {', '.join(supported_solvents.keys())}"
        )

    return supported_solvents[solvent_name]


def validate_input(smiles, solvent_name, supported_solvents):
    """
    Validate prediction inputs.

    Args:
        smiles (str): molecule SMILES
        solvent_name (str): solvent name
        supported_solvents (dict): dictionary of supported solvents

    Returns:
        tuple: (canonical_smiles, solvent_smiles)

    Raises:
        InvalidSMILESError: if SMILES is invalid
        ValidationError: if solvent is invalid
    """
    # Validate SMILES
    validate_smiles(smiles)

    # Canonicalize SMILES
    mol = Chem.MolFromSmiles(smiles)
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    # Validate solvent
    solvent_smiles = validate_solvent(solvent_name, supported_solvents)

    return canonical_smiles, solvent_smiles
