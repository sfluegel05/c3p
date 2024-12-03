"""
Classifies: CHEBI:26469 quaternary nitrogen compound
"""
from rdkit import Chem

def is_quaternary_nitrogen_compound(smiles: str):
    """
    Determines if a molecule is a quaternary nitrogen compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary nitrogen compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    quaternary_nitrogen_found = False

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check if the atom is nitrogen
            if atom.GetFormalCharge() == 1:  # Check if the nitrogen has a positive charge
                neighbors = atom.GetNeighbors()
                if len(neighbors) == 4:  # Check if the nitrogen is bonded to four other atoms
                    quaternary_nitrogen_found = True
                    break

    if quaternary_nitrogen_found:
        return True, "Quaternary nitrogen found"
    else:
        return False, "No quaternary nitrogen found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26469',
                          'name': 'quaternary nitrogen compound',
                          'definition': 'A nitrogen molecular entity that is '
                                        'electronically neutral but which '
                                        'contains a quaternary nitrogen.',
                          'parents': ['CHEBI:35352']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[15:22:41] Explicit valence for atom # 1 N, 5, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 200,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 8,
    'precision': 1.0,
    'recall': 0.9615384615384616,
    'f1': 0.9803921568627451,
    'accuracy': None}