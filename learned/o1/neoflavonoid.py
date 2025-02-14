"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for coumarin core with aryl substituent at position 4
    # Coumarin core (1-benzopyran-2-one) with atom mapping
    # Atom indices:
    # [O:1]=[C:2]1[C:3]=[C:4][C:5]=[C:6][C:7]=[C:8][O:9]1
    # Aryl substituent attached to [C:3]
    neoflavonoid_smarts = "[O]=C1C=CC=CC(O1)[C]-[a]"
    neoflavonoid_pattern = Chem.MolFromSmarts(neoflavonoid_smarts)

    if neoflavonoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the neoflavonoid pattern in the molecule
    match = mol.HasSubstructMatch(neoflavonoid_pattern)
    if match:
        return True, "Molecule is a neoflavonoid with aryl substituent at position 4 of the coumarin core"
    else:
        return False, "Neoflavonoid core with aryl substituent at position 4 not found"