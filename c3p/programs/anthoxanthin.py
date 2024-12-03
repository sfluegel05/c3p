"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a flavonoid
    flavonoid_basic_structure = Chem.MolFromSmarts("c1cc2oc(c3ccc(O)cc3)c(=O)c2c(O)c1")
    if not mol.HasSubstructMatch(flavonoid_basic_structure):
        return False, "Molecule does not match the basic flavonoid structure"

    # Check for additional hydroxyl groups
    hydroxyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]"))
    if len(hydroxyl_groups) < 2:
        return False, "Molecule does not have sufficient hydroxyl groups"

    # Check for methoxy groups
    methoxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("CO"))
    if len(methoxy_groups) == 0:
        return False, "Molecule does not have any methoxy groups"

    # Check for glycosylation (sugar moieties)
    sugar_smarts = Chem.MolFromSmarts("[C@H]1(O[C@H](O[C@H](O[C@H]1O)CO)CO)O")
    if not mol.HasSubstructMatch(sugar_smarts):
        return False, "Molecule does not have glycosylation"

    return True, "Molecule matches the anthoxanthin structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:192499',
                          'name': 'anthoxanthin',
                          'definition': 'are a type of flavonoid pigments in '
                                        'plants. Anthoxanthins are '
                                        'water-soluble pigments which range in '
                                        'color from white or colorless to a '
                                        'creamy to yellow, often on petals of '
                                        'flowers.',
                          'parents': ['CHEBI:47916']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 66,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}