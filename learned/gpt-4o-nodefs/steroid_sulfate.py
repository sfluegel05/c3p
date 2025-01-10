"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.

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

    # Look for sulfate groups (OS(=O)(=O)O)
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-] | OS(O)(=O)=O")
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate groups found"

    # Look for steroid backbone (4 fused rings structure)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C4CCC(C4)C3CC2C1")  # This pattern represents typical steroid-like structure
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    return True, "Contains steroid backbone with sulfate group(s)"