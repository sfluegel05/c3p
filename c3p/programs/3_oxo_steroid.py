"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid (Any oxo steroid where an oxo substituent is located at position 3).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a steroid backbone
    steroid_backbone = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCC4')  # Simplified steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for the presence of an oxo group at position 3
    oxo_group = Chem.MolFromSmarts('C=O')
    if not mol.HasSubstructMatch(oxo_group):
        return False, "No oxo group found"

    # Find the oxo group and check its position
    matches = mol.GetSubstructMatches(oxo_group)
    for match in matches:
        carbonyl_carbon = match[0]
        neighbors = mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon
                if neighbor.GetIdx() == 2:  # Check if the position is 3 (0-indexed)
                    return True, "3-oxo steroid found"

    return False, "No oxo group at position 3 found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47788',
                          'name': '3-oxo steroid',
                          'definition': 'Any oxo steroid where an oxo '
                                        'substituent is located at position 3.',
                          'parents': ['CHEBI:35789', 'CHEBI:3992']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 25-26: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}