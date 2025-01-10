"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

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

    # Look for sulfate groups
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[OH]")
    charged_sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")
    
    sulfate_match = mol.HasSubstructMatch(sulfate_pattern)
    charged_sulfate_match = mol.HasSubstructMatch(charged_sulfate_pattern)
    
    if not (sulfate_match or charged_sulfate_match):
        return False, "No sulfate groups found"

    # Look for steroid backbone (fused ring system)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C4CCC(C4)C3CC2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    return True, "Contains steroid backbone with sulfate group(s)"