"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phytosphingosine backbone
    phytosphingosine_pattern = Chem.MolFromSmarts('O[C@H](CO)C[C@H](O)C[C@H](O)CCCCCCCCCCCCCC')
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Check for the presence of an acyl group attached to the nitrogen
    acyl_group_pattern = Chem.MolFromSmarts('NC(=O)C')
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No acyl group attached to the nitrogen found"

    # Check for the presence of a glycosyl group attached to oxygen (optional)
    glycosyl_group_pattern = Chem.MolFromSmarts('O[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O')
    if mol.HasSubstructMatch(glycosyl_group_pattern):
        return True, "N-acylphytosphingosine with glycosyl group"

    return True, "N-acylphytosphingosine without glycosyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:31998',
                          'name': 'N-acylphytosphingosine',
                          'definition': 'A ceramide that is phytosphingosine '
                                        'having a fatty acyl group attached to '
                                        'the nitrogen.',
                          'parents': ['CHEBI:139051']},
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