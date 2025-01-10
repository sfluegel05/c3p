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

    # Define a general SMARTS pattern for the butenolide ring
    # This pattern matches a five-membered ring with one oxygen atom,
    # one carbonyl group (C=O), and one double bond, allowing substitutions
    butenolide_smarts = '[#8]=C1[#6][#6]=[#6][#8]1'  # Generalized 2-furanone skeleton
    butenolide_pattern = Chem.MolFromSmarts(butenolide_smarts)

    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Contains butenolide ring (2-furanone skeleton with substitutions)"
    else:
        return False, "No butenolide ring found"