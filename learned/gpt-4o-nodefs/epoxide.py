"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    Epoxides are three-membered cyclic ethers, typically represented by the pattern C1OC1.

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

    # Define the epoxide pattern (3-membered ring with an oxygen atom)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains a 3-membered cyclic ether (epoxide)"
    else:
        return False, "No epoxide structure detected"