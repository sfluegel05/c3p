"""
Classifies: CHEBI:26961 thiophenes
"""
from rdkit import Chem

def is_thiophenes(smiles: str):
    """
    Determines if a molecule contains at least one thiophene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains at least one thiophene ring, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the thiophene substructure
    thiophene_smarts = 'c1ccsc1'
    thiophene = Chem.MolFromSmarts(thiophene_smarts)

    # Check if the molecule contains the thiophene substructure
    if mol.HasSubstructMatch(thiophene):
        return True, "Contains at least one thiophene ring"
    else:
        return False, "Does not contain a thiophene ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26961',
                          'name': 'thiophenes',
                          'definition': 'Compounds containing at least one '
                                        'thiophene ring.',
                          'parents': ['CHEBI:38106']},
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
    'num_true_positives': 33,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}