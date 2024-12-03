"""
Classifies: CHEBI:76170 olefinic phospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_olefinic_phospholipid(smiles: str):
    """
    Determines if a molecule is an olefinic phospholipid (phospholipid containing an olefinic function).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic phospholipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phospholipid by presence of phosphoric acid ester group
    phosphoric_acid_ester = Chem.MolFromSmarts('COP(O)(=O)O')
    if not mol.HasSubstructMatch(phosphoric_acid_ester):
        return False, "No phosphoric acid ester group found"

    # Check for olefinic function (C=C double bond)
    olefinic_function = Chem.MolFromSmarts('C=C')
    if not mol.HasSubstructMatch(olefinic_function):
        return False, "No olefinic function (C=C double bond) found"

    # Check for the presence of a glycerol backbone
    glycerol_backbone = Chem.MolFromSmarts('[C@H](CO)O')
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    return True, "Molecule is an olefinic phospholipid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76170',
                          'name': 'olefinic phospholipid',
                          'definition': 'Any phospholipid that contains an '
                                        'olefinic function.',
                          'parents': ['CHEBI:16247', 'CHEBI:78840']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 2,
    'num_true_negatives': 13,
    'num_false_negatives': 0,
    'precision': 0.8823529411764706,
    'recall': 1.0,
    'f1': 0.9375,
    'accuracy': None}