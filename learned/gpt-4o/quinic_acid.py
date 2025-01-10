"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or a derivative based on its SMILES string.
    Quinic acid has a cyclitol skeletal structure with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as quinic acid or derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclohexane 6-membered ring
    cyclohexane_pattern = Chem.MolFromSmarts("C1CCCCC1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane ring found"

    # Look for at least three hydroxyl groups on the cyclohexane
    hydroxyl_on_ring_pattern = Chem.MolFromSmarts("[C;R1][OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_on_ring_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Only {len(hydroxyl_matches)} hydroxyl groups found on cyclohexane, need at least 3"

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Optional: check for ester linkage to aromatic acids (e.g., caffeoyl or coumaroyl)
    aromatic_ester_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)O[C;R1]")
    if mol.HasSubstructMatch(aromatic_ester_pattern):
        return True, "Quinic acid derivative with aromatic ester linkage found"

    return True, "Quinic acid detected with basic structural features"