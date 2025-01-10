"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:32994 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a cyclic ether in which the oxygen atom forms part of a 3-membered ring.

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

    # Define the epoxide pattern: a 3-membered ring with one oxygen and two carbons
    epoxide_pattern = Chem.MolFromSmarts("[O;R3][C;R3][C;R3]")
    
    # Check if the molecule contains the epoxide pattern
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains a 3-membered ring with an oxygen atom (epoxide)"
    else:
        return False, "No 3-membered ring with an oxygen atom found"