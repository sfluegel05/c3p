"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: 2-hydroxy fatty acid (CHEBI:78283)
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES.
    A 2-hydroxy fatty acid has a carboxylic acid group with a hydroxyl on the adjacent (alpha/2-position) carbon,
    and a sufficiently long carbon chain typical of fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if criteria are met, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for at least one carboxylic acid group
    carboxyl_pattern = MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group"

    # Check for hydroxyl on alpha carbon (adjacent to carboxylic acid's carbonyl)
    alpha_hydroxy_pattern = MolFromSmarts("[CX3](=O)-[CH1]-[OH]")  # Alpha carbon with hydroxyl
    alpha_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    if not alpha_matches:
        return False, "No hydroxyl on alpha carbon"

    # Verify sufficient chain length (total carbons >=5 to exclude small acids like glycolic)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, f"Insufficient carbons ({c_count}) for fatty acid"

    return True, "Carboxylic acid with hydroxyl on alpha carbon and sufficient chain length"