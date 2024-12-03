"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import Descriptors


def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin (a saponin derived from a hydroxysteroid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a steroid backbone
    steroid_backbone = Chem.MolFromSmarts('C1CCC2C3CCC4CC5CCCCC5C4C3C2C1')
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for hydroxyl groups attached to the steroid backbone
    hydroxyl_group = Chem.MolFromSmarts('[C,c][OH]')
    if not mol.HasSubstructMatch(hydroxyl_group):
        return False, "No hydroxyl groups found attached to the steroid backbone"

    # Check for glycosidic linkages (sugar moieties)
    glycosidic_linkage = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)[C@H](O)[C@H]1O')
    if not mol.HasSubstructMatch(glycosidic_linkage):
        return False, "No glycosidic linkages found"

    return True, "Molecule is a steroid saponin"

# Example usage
smiles = 'C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H](O)C4=C([C@H]3CC[C@]12C)C(=O)C[C@H](O)C4)[C@H]1CC(CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)=C(C)C(=O)O1'
print(is_steroid_saponin(smiles))  # Expected output: (True, "Molecule is a steroid saponin")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61655',
                          'name': 'steroid saponin',
                          'definition': 'Any saponin derived from a '
                                        'hydroxysteroid.',
                          'parents': ['CHEBI:26605', 'CHEBI:35341']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'No steroid backbone found')\n",
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 36,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}