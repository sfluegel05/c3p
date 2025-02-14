"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is a nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitrile group (C#N)
    nitrile_pattern = Chem.MolFromSmarts("C#N")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile group found"
    
    # Check for aromatic rings
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Aromatic ring detected, not an aliphatic nitrile."

    return True, "Contains a nitrile group and no aromatic rings."