"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the butenolide pattern (non-aromatic form)
    butenolide_smarts = 'O=C1OC=CC1'  # 2-furanone skeleton
    butenolide_pattern = Chem.MolFromSmarts(butenolide_smarts)

    # Define the aromatic butenolide pattern
    aromatic_butenolide_smarts = 'O=c1occc1'  # Aromatic form
    aromatic_butenolide_pattern = Chem.MolFromSmarts(aromatic_butenolide_smarts)

    # Check for the butenolide ring (non-aromatic)
    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Contains butenolide ring (2-furanone skeleton)"
    # Check for the aromatic butenolide ring
    elif mol.HasSubstructMatch(aromatic_butenolide_pattern):
        return True, "Contains aromatic butenolide ring (2-furanone skeleton)"
    else:
        return False, "No butenolide ring found"