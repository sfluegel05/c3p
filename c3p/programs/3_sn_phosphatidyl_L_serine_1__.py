"""
Classifies: CHEBI:57262 3-sn-phosphatidyl-L-serine(1-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_3_sn_phosphatidyl_L_serine_1__(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine(1-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine(1-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphoserine group
    phosphoserine_pattern = Chem.MolFromSmarts('[NH3+][C@@H](COP([O-])(=O)OC[C@@H](COC([*])=O)OC([*])=O)C([O-])=O')
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # Check for the presence of the glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts('OCC(O)CO')
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for the presence of the diacyl groups
    diacyl_pattern = Chem.MolFromSmarts('C(=O)OCC(=O)O')
    if not mol.HasSubstructMatch(diacyl_pattern):
        return False, "No diacyl groups found"

    return True, "Molecule is a 3-sn-phosphatidyl-L-serine(1-)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:57262',
                          'name': '3-sn-phosphatidyl-L-serine(1-)',
                          'definition': 'A singly-charged anionic phospholipid '
                                        'that is the conjugate base of '
                                        '1,2-diacyl-sn-glycero-3-phosphoserine, '
                                        'in which the carboxy and phosphate '
                                        'groups are anionic and the amino '
                                        'group is cationic; major species at '
                                        'pH 7.3.',
                          'parents': ['CHEBI:58944', 'CHEBI:62643']},
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
    'num_true_negatives': 11,
    'num_false_negatives': 11,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}