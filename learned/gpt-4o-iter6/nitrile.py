"""
Classifies: CHEBI:18379 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a cyano group (-C#N), which is the defining feature.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyano group pattern (C#N)
    cyano_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    if not mol.HasSubstructMatch(cyano_pattern):
        return False, "No cyano group (C#N) found"
    
    return True, "Contains cyano group (C#N), identified as a nitrile"