"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid typically contains a backbone with a hydroxyl group on the
    second carbon and two carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the 2-hydroxy group (OH on the second carbon)
    hydroxy_pattern = Chem.MolFromSmarts("[CH2X4;R0][CHX4;R0][OX2H]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 2-hydroxy group found"

    # Look for two carboxylic acid groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 2:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 2"

    return True, "Contains 2-hydroxy group and two carboxylic acid groups"