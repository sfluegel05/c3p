"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:32987 epoxide
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

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an epoxide
    epoxide_pattern = Chem.MolFromSmarts('[OX2r3][CX4r3][CX4r3]')

    # Check if the molecule matches the epoxide pattern
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains an epoxide ring (3-membered ring with one oxygen atom)"
    else:
        return False, "Does not contain an epoxide ring"