"""
Classifies: CHEBI:36700 phosphocholines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphocholines(smiles: str):
    """
    Determines if a molecule is a phosphocholine (any compound having phosphocholine as part of its structure).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphocholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the phosphocholine substructure
    phosphocholine_smarts = 'COP(OCC[N+](C)(C)C)(=O)[O-]'
    phosphocholine = Chem.MolFromSmarts(phosphocholine_smarts)

    if phosphocholine is None:
        return None, "Failed to create phosphocholine substructure"

    # Check if the molecule contains the phosphocholine substructure
    if mol.HasSubstructMatch(phosphocholine):
        return True, "Molecule contains phosphocholine substructure"
    else:
        return False, "Molecule does not contain phosphocholine substructure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36700',
                          'name': 'phosphocholines',
                          'definition': 'Any compound having phosphocholine as '
                                        'part of its structure.',
                          'parents': [   'CHEBI:23213',
                                         'CHEBI:23217',
                                         'CHEBI:25703',
                                         'CHEBI:37734']},
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
    'num_true_positives': 212,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.986046511627907,
    'f1': 0.9929742388758782,
    'accuracy': None}