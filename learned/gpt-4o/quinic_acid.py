"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or one of its derivatives based on its SMILES string.
    Quinic acid is characterized as a cyclitol carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as quinic acid or a derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pattern for a cyclohexane ring with at least 3 hydroxyl groups
    hydroxylated_cyclohexane_pattern = Chem.MolFromSmarts("[C@H]1[C@H][C@H][C@H][C@H][C@H]1")
    hydroxyl_groups_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([OH])[C@H]([OH])[C@H]([OH])[C@H][C@H]1")
    if not mol.HasSubstructMatch(hydroxylated_cyclohexane_pattern):
        return False, "No cyclohexane ring with correct stereochemistry found"

    # Ensure there are at least 3 hydroxyl groups on the cyclohexane ring
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_groups_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "Less than three hydroxyl groups on cyclohexane ring with required stereochemistry"

    # Look for a carboxylic acid group potentially linked to the cyclohexane
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Detect ester linkages specifically for quinic acid derivatives
    ester_pattern = Chem.MolFromSmarts("c1ccc(O[C@H]2[C@H](O)[C@@H](O)[C@@H]([C@@H]2O)C(=O)O)cc1")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Quinic acid derivative with an aromatic ester linkage found"

    return True, "Quinic acid detected with basic structural features"