"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is based on 1-benzopyran (chromene) with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for flavonoid cores
    flavonoid_smarts_list = [
        # Flavone core
        '[O]=C1C=CC(=CC1)-c1ccccc1',
        # Flavonol core
        '[O]=C1C=CC(=CC1O)-c1ccccc1',
        # Flavanone core
        'O=C1CC=C(C=C1)-c1ccccc1',
        # Flavanol (catechin) core
        'O1CC=C(C=C1)-c1cccc(O)c1',
        # Isoflavone core
        '[O]=C1C=CC(=CC1-c1ccccc1)',
        # Anthocyanidin core
        '[O+]C1=CC=C(O)C=C1-c1cccc(O)c1',
    ]

    # Check if molecule matches any of the flavonoid cores
    for flavonoid_smarts in flavonoid_smarts_list:
        pattern = Chem.MolFromSmarts(flavonoid_smarts)
        if pattern is None:
            continue  # Skip invalid patterns
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule matches flavonoid core structure"
    
    return False, "Molecule does not match flavonoid core structure"