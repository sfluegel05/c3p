"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide
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

    # Define the butenolide SMARTS pattern (2-furanone skeleton)
    butenolide_smarts = "O=C1C=CCO1"
    pattern = Chem.MolFromSmarts(butenolide_smarts)

    # Check if molecule contains the butenolide core
    if mol.HasSubstructMatch(pattern):
        return True, "Contains butenolide ring (2-furanone skeleton)"
    else:
        return False, "Does not contain butenolide ring (2-furanone skeleton)"