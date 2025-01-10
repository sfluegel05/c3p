"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as any aldehyde with the carbonyl group attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the aldehyde group (C=O attached to one hydrogen and another carbon)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"

    # Look for specific aromatic ring patterns attached to the aldehyde's carbon
    aromatic_patterns = [
        Chem.MolFromSmarts("[c]"),  # aromatic carbon
        Chem.MolFromSmarts("n1ccccc1"),  # pyridine-like
        Chem.MolFromSmarts("[cH][c][c][c][c][c]"),  # benzene-like
        Chem.MolFromSmarts("n2c[nH]nc2"),  # pyrrole-like with nitrogen consideration
    ]
    
    # Iterate through matched aldehyde carbons to check for attachment to an aromatic system
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    for match in matches:
        aldehyde_carbon = match[0]
        for neighbor in mol.GetAtomWithIdx(aldehyde_carbon).GetNeighbors():
            # Check if the neighbor is part of any of the defined aromatic patterns
            for aromatic_pattern in aromatic_patterns:
                if mol.HasSubstructMatch(aromatic_pattern, useChirality=False):
                    return True, "Aldehyde group attached to an aromatic moiety"

    return False, "Aldehyde group not attached to aromatic moiety"