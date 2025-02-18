"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion has a central nitrogen atom with four organic substituents and a positive charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for nitrogen with four carbon substituents and positive charge
    quaternary_n_pattern = Chem.MolFromSmarts("[N+]([#6])([#6])([#6])[#6]")
    # Pattern to detect carboxylate groups (common in false positives like carnitine)
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")

    # Check for quaternary nitrogen and absence of carboxylate groups
    if mol.HasSubstructMatch(quaternary_n_pattern):
        if not mol.HasSubstructMatch(carboxylate_pattern):
            return True, "Contains a nitrogen atom with four organic substituents and a positive charge"
        else:
            return False, "Quaternary nitrogen present but molecule contains carboxylate group (common false positive)"
    
    return False, "No quaternary ammonium group found"