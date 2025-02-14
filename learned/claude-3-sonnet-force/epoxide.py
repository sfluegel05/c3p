"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:33569 epoxide
An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.

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
    
    # Look for epoxide ring pattern (O linked to 2 carbons with a single bond each)
    epoxide_pattern = Chem.MolFromSmarts("[O;r3]")
    
    # Check if the molecule contains at least one match for the epoxide pattern
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Molecule contains an epoxide ring (a 3-membered cyclic ether)"
    else:
        return False, "No epoxide ring found in the molecule"