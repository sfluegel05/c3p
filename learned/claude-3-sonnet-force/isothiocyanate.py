"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: CHEBI:50974 isothiocyanate
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate is an organosulfur compound with the general formula R-N=C=S.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isothiocyanate functional group pattern (N=C=S)
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
    if not mol.HasSubstructMatch(isothiocyanate_pattern):
        return False, "No isothiocyanate functional group (N=C=S) found"
    
    # Check for carbon-nitrogen bond (R-N=C=S)
    has_c_n_bond = any(bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 7
                       for bond in mol.GetBonds())
    if not has_c_n_bond:
        return False, "No carbon-nitrogen bond found (R-N=C=S)"
    
    return True, "Contains isothiocyanate functional group (R-N=C=S)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:50974',
        'name': 'isothiocyanate',
        'definition': 'An organosulfur compound with the general formula R-N=C=S.',
        'parents': ['CHEBI:24432']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 175,
    'num_false_positives': 1,
    'num_true_negatives': 182414,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.9942857142857143,
    'recall': 1.0,
    'f1': 0.9971098265895953,
    'accuracy': 0.9999994547617611
}