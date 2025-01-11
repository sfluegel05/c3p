"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide substructure (C1OC1 pattern)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains the 3-membered ether ring indicating an epoxide"

    return False, "Does not contain the 3-membered ether ring characteristic of an epoxide"