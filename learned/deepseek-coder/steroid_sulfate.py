"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: CHEBI:37943 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    steroid_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@](CC1)(CC[C@@H]4[C@@H]3CC2)C)C)")  # General steroid backbone
    sulfate_pattern = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX1H0-]")  # Sulfate group

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for sulfate group
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group found"

    # Ensure the sulfate group is attached to the steroid backbone
    steroid_atoms = set(mol.GetSubstructMatch(steroid_pattern))
    sulfate_atoms = set(mol.GetSubstructMatch(sulfate_pattern))
    if not steroid_atoms.intersection(sulfate_atoms):
        return False, "Sulfate group not attached to steroid backbone"

    return True, "Contains steroid backbone with attached sulfate group"