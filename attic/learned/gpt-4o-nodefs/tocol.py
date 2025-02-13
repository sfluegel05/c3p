"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tocol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for chroman or chromene moiety - ring structure with an oxygen and aromatic characteristics
    chroman_pattern = Chem.MolFromSmarts("O1C(C=C)=C(CC(C)C)C(CCC(C)(C)C)C1")
    if not mol.HasSubstructMatch(chroman_pattern):
        return False, "No chroman or chromene moiety"

    # Check if there's a long hydrocarbon tail
    hydrocarbon_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(hydrocarbon_pattern):
        return False, "No long hydrocarbon chain found"

    # Check for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl groups"

    return True, "Contains chroman/chromene moiety with hydrocarbon tail and hydroxyl groups"