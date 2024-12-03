"""
Classifies: CHEBI:50402 androstanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_androstanoid(smiles: str):
    """
    Determines if a molecule is an androstanoid (based on an androstane skeleton and its derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an androstanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the androstane skeleton SMARTS pattern
    androstane_skeleton = Chem.MolFromSmarts('C1CCC2C3CCC4CCCC[C@]4(C)C3CCC12')
    
    if mol.HasSubstructMatch(androstane_skeleton):
        return True, "Molecule contains the androstane skeleton"
    else:
        return False, "Molecule does not contain the androstane skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50402',
                          'name': 'androstanoid',
                          'definition': 'Any steroid based on an  androstane '
                                        'skeleton and its derivatives.',
                          'parents': ['CHEBI:35341']},
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
    'num_true_positives': 7,
    'num_false_positives': 2,
    'num_true_negatives': 9,
    'num_false_negatives': 4,
    'precision': 0.7777777777777778,
    'recall': 0.6363636363636364,
    'f1': 0.7000000000000001,
    'accuracy': None}