"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as a monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    The ring might also have other non-carbon atoms such as sulfur or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for common azoles
    azole_patterns = [
        Chem.MolFromSmarts("n1cc[n,o,s,c]c1"),  # Generic azole, any heteroatom
        Chem.MolFromSmarts("c1n[c,n,o][c,n,o,s]c1"),  # Variant with carbon/nitrogen positions
    ]

    # Check each azole pattern
    for pattern in azole_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains azole ring (5-membered ring with nitrogen and other heteroatoms)"

    return False, "Does not contain an azole ring"