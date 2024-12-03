"""
Classifies: CHEBI:33854 aromatic alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_aromatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an aromatic ring
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('a')):
        return False, "No aromatic ring found"

    # Check for the presence of an alcohol group (-OH)
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[CX4][OH]')):
        return False, "No alcohol group found"

    # Check if the alcohol group is attached to a carbon which is bonded to an aromatic ring
    alcohol_smarts = Chem.MolFromSmarts('[CX4][OH]')
    matches = mol.GetSubstructMatches(alcohol_smarts)

    for match in matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Aromatic alcohol found"

    return False, "Alcohol group is not attached to a carbon bonded to an aromatic ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33854',
                          'name': 'aromatic alcohol',
                          'definition': 'Any alcohol in which the alcoholic '
                                        'hydroxy group is attached to a carbon '
                                        'which is itself bonded to an aromatic '
                                        'ring.',
                          'parents': ['CHEBI:30879', 'CHEBI:33659']},
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
    'num_true_positives': 16,
    'num_false_positives': 2,
    'num_true_negatives': 15,
    'num_false_negatives': 1,
    'precision': 0.8888888888888888,
    'recall': 0.9411764705882353,
    'f1': 0.9142857142857143,
    'accuracy': None}