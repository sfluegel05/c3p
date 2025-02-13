"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound, i.e.,
    it contains at least one hydroxyl group (-OH) attached to a non-aromatic carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found"

    # Look for aliphatic alcohol hydroxyl groups
    aliphatic_alcohol_pattern = Chem.MolFromSmarts("[C;!a][OX2H]")
    aliphatic_alcohol_matches = mol.GetSubstructMatches(aliphatic_alcohol_pattern)
    if len(aliphatic_alcohol_matches) == 0:
        return False, "No aliphatic hydroxyl groups found (only aromatic hydroxyl groups)"
    else:
        return True, "Contains aliphatic hydroxyl group(s)"

__metadata__ = {
    'chemical_class': {
        'name': 'aliphatic alcohol',
        'definition': 'An alcohol derived from an aliphatic compound.',
    }
}