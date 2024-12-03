"""
Classifies: CHEBI:62941 glycosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glycosylceramide(smiles: str):
    """
    Determines if a molecule is a glycosylceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosylceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a ceramide group
    ceramide_pattern = Chem.MolFromSmarts('[C@H](O)[C@H](NC(=O)*)CO')
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide group found"

    # Check for the presence of a glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H]([C@H](O)[C@H](O)[C@H](O1)CO)')
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found"

    # Check for the presence of a monosaccharide
    monosaccharide_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H]([C@H](O)[C@H](O)[C@H](O1)CO)')
    if not mol.HasSubstructMatch(monosaccharide_pattern):
        return False, "No monosaccharide found"

    return True, "Molecule is a glycosylceramide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62941',
                          'name': 'glycosylceramide',
                          'definition': 'A ceramide compound formed by the '
                                        'replacement of the glycosidic hydroxy '
                                        'group of a cyclic form of a '
                                        'monosaccharide (or derivative) by a '
                                        'ceramide group.',
                          'parents': ['CHEBI:17761', 'CHEBI:24402']},
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
    'num_true_positives': 27,
    'num_false_positives': 17,
    'num_true_negatives': 3,
    'num_false_negatives': 0,
    'precision': 0.6136363636363636,
    'recall': 1.0,
    'f1': 0.7605633802816901,
    'accuracy': None}