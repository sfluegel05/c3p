"""
Classifies: CHEBI:51751 enamide
"""
from rdkit import Chem

def is_enamide(smiles: str):
    """
    Determines if a molecule is an enamide (alpha,beta-unsaturated carboxylic acid amide).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enamide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for enamide
    enamide_smarts = 'C=CC(=O)NC'
    
    # Check if the molecule matches the enamide pattern
    enamide_pattern = Chem.MolFromSmarts(enamide_smarts)
    if mol.HasSubstructMatch(enamide_pattern):
        return True, "Molecule is an enamide"
    else:
        return False, "Molecule is not an enamide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51751',
                          'name': 'enamide',
                          'definition': 'An alpha,beta-unsaturated carboxylic '
                                        'acid amide of general formula '
                                        'R(1)R(2)C=CR(3)-C(=O)NR(4)R(5) in '
                                        'which the amide C=O function is '
                                        'conjugated to a C=C double bond at '
                                        'the alpha,beta position.',
                          'parents': ['CHEBI:51750', 'CHEBI:78840']},
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
    'num_true_positives': 14,
    'num_false_positives': 1,
    'num_true_negatives': 16,
    'num_false_negatives': 3,
    'precision': 0.9333333333333333,
    'recall': 0.8235294117647058,
    'f1': 0.8749999999999999,
    'accuracy': None}