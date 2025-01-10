"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile group (C#N) connected to an aliphatic (non-aromatic)
    carbon chain or environment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrile group pattern (-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    
    if not nitrile_matches:
        return False, "No nitrile group found"

    # A SMARTS pattern ensuring connection to non-aromatic systems
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CX4,CX3][CX2]#N")
    
    # Enhanced check focused on connection to non-aromatic systems directly:
    for match in nitrile_matches:
        for path in mol.GetSubstructMatches(aliphatic_chain_pattern):
            if match == path[1:]:  # If nitrile is part of such a pattern
                return True, "Contains an aliphatic nitrile group"

    return False, "Nitrile group is not connected to an aliphatic chain"