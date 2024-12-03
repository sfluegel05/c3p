"""
Classifies: CHEBI:25036 lignan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_lignan(smiles: str):
    """
    Determines if a molecule is a lignan.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lignan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of dibenzylbutane skeleton
    try:
        # SMARTS pattern for dibenzylbutane skeleton
        pattern = Chem.MolFromSmarts('C1CC(C2=CC=CC=C2)C(C3=CC=CC=C3)C1')
        
        if mol.HasSubstructMatch(pattern):
            return True, "Contains dibenzylbutane skeleton"
        else:
            return False, "Does not contain dibenzylbutane skeleton"
    except:
        return None, None

    return False, "Unknown reason"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25036',
                          'name': 'lignan',
                          'definition': 'Any phenylpropanoid  derived from '
                                        'phenylalanine via dimerization of '
                                        'substituted cinnamic alcohols, known '
                                        'as monolignols, to a dibenzylbutane '
                                        'skeleton. Note that while individual '
                                        'members of the class have names '
                                        'ending ...lignane, ...lignene, '
                                        '...lignadiene, etc., the class names '
                                        'lignan, neolignan, etc., do not end '
                                        'with an "e".',
                          'parents': ['CHEBI:26004']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 38,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}