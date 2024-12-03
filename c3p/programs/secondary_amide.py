"""
Classifies: CHEBI:33257 secondary amide
"""
from rdkit import Chem

def is_secondary_amide(smiles: str):
    """
    Determines if a molecule is a secondary amide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the secondary amide functional group
    amide_patterns = [
        Chem.MolFromSmarts('C(=O)N(C)C(=O)'),
        Chem.MolFromSmarts('C(=O)N(C)'),
        Chem.MolFromSmarts('N(C)C(=O)')
    ]
    
    for pattern in amide_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a secondary amide group"
    
    return False, "No secondary amide group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33257',
                          'name': 'secondary amide',
                          'definition': 'A derivative of two oxoacids '
                                        'RkE(=O)l(OH)m (l =/= 0) in which two '
                                        'acyl groups are attached to the amino '
                                        'or substituted amino group.',
                          'parents': ['CHEBI:32988']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 23,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 2,
    'precision': 0.7419354838709677,
    'recall': 0.92,
    'f1': 0.8214285714285714,
    'accuracy': None}