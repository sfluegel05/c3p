"""
Classifies: CHEBI:68489 alkyl,acyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alkyl_acyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is an alkyl, acyl-sn-glycero-3-phosphocholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkyl, acyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts('COP([O-])(=O)OCC[N+](C)(C)C')
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for the glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts('C(CO)O')
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for the presence of at least one alkyl and one acyl group
    alkyl_pattern = Chem.MolFromSmarts('C(C)C')
    acyl_pattern = Chem.MolFromSmarts('C(=O)C')

    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)

    if len(alkyl_matches) < 1:
        return False, "No alkyl group found"
    if len(acyl_matches) < 1:
        return False, "No acyl group found"

    return True, "Molecule is an alkyl, acyl-sn-glycero-3-phosphocholine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:68489',
                          'name': 'alkyl,acyl-sn-glycero-3-phosphocholine',
                          'definition': 'A glycerophosphocholine that is '
                                        'sn-glycero-3-phosphocholine in which '
                                        'positions 1 and 2 are substituted by '
                                        'unspecified acyl and alkyl groups, '
                                        'and in which the positions of the '
                                        'acyl and alkyl groups are also '
                                        'unspecified.',
                          'parents': ['CHEBI:36313', 'CHEBI:78189']},
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
    'num_true_positives': 9,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 13,
    'precision': 0.6,
    'recall': 0.4090909090909091,
    'f1': 0.4864864864864865,
    'accuracy': None}