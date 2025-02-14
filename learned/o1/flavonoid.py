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
    A flavonoid is based on a 1-benzopyran (chromene) with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the flavonoid core:
    # - A benzopyran ring system (fusion of benzene and pyran rings)
    # - An aryl group attached at position 2 of the pyran ring
    flavonoid_smarts = """
    [$([cH]1[cH][cH][cH][cH][cH]1)]               # Benzene ring (A ring)
    [C]                                            # Carbon connecting A and C rings
    2                                             # Ring closure to form C ring
    [O][cH][cH][cH]2                               # Pyran ring (C ring) with oxygen and closure
    [c]                                            # Aromatic carbon at position 2 of pyran ring
    [c]:[c]:[c]:[c]:[c]                            # B ring (aryl substituent)
    """

    # Remove whitespace and newlines from SMARTS pattern
    flavonoid_smarts = ''.join(flavonoid_smarts.strip().split())

    # Create SMARTS pattern
    pattern = Chem.MolFromSmarts(flavonoid_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for flavonoid core match
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule matches flavonoid core structure"

    return False, "Molecule does not match flavonoid core structure"