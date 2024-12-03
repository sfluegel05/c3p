"""
Classifies: CHEBI:38102 triazines
"""
from rdkit import Chem

def is_triazines(smiles: str):
    """
    Determines if a molecule is a triazine (compounds based on a triazine skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triazine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS patterns for different types of triazine rings
    triazine_patterns = [
        Chem.MolFromSmarts('c1ncncn1'),  # 1,3,5-triazine
        Chem.MolFromSmarts('c1ncnnc1'),  # 1,2,4-triazine
        Chem.MolFromSmarts('c1ncn[nH]c1'),  # 1,2,3-triazine
        Chem.MolFromSmarts('n1cncnc1')  # 1,2,3,5-tetrazine
    ]
    
    for pattern in triazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a triazine ring"
    
    return False, "Molecule does not contain a triazine ring"

# Example usage:
# print(is_triazines("Nc1ncn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1"))
# Output: (True, "Molecule contains a triazine ring")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38102',
                          'name': 'triazines',
                          'definition': 'Compounds based on a triazine '
                                        'skeleton.',
                          'parents': [   'CHEBI:25693',
                                         'CHEBI:38101',
                                         'CHEBI:50893']},
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
    'num_true_positives': 17,
    'num_false_positives': 0,
    'num_true_negatives': 18,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9444444444444444,
    'f1': 0.9714285714285714,
    'accuracy': None}