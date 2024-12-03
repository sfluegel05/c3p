"""
Classifies: CHEBI:33250 atom
"""
from rdkit import Chem

def is_atom(smiles: str):
    """
    Determines if a molecule is an atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an atom, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    if mol.GetNumAtoms() == 1:
        return True, "The SMILES string represents a single atom"
    else:
        return False, "The SMILES string does not represent a single atom"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33250',
                          'name': 'atom',
                          'definition': 'A chemical entity constituting the '
                                        'smallest component of an element '
                                        'having the chemical properties of the '
                                        'element.',
                          'parents': ['CHEBI:24431']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[22:59:07] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 23,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}