"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: phytosterols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone SMARTS pattern (tetracyclic ring system)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C1CCC3C2CC=C4C3CCCC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for hydroxyl group at position 3 (3-beta OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[C@H](O)[C@@H]1CCCC2=C1CCCC2')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group at position 3"

    # Check that the molecule is not cholesterol itself
    cholesterol_smiles = "C[C@H](CCC(=C)C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)O)C"
    cholesterol_mol = Chem.MolFromSmiles(cholesterol_smiles)
    if mol.HasSubstructMatch(cholesterol_mol):
        return False, "Molecule is cholesterol, not a phytosterol"

    # If all checks passed, classify as phytosterol
    return True, "Molecule is a phytosterol with steroid backbone and variations in side chain/double bonds"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'phytosterols',
        'definition': 'Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here if needed
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata fields can be added here
}