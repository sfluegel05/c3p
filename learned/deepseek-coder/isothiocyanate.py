"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

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

    # Look for the isothiocyanate group (-N=C=S)
    isothiocyanate_pattern = Chem.MolFromSmarts("[N]=[C]=[S]")
    if not mol.HasSubstructMatch(isothiocyanate_pattern):
        return False, "No isothiocyanate group (-N=C=S) found"

    # Check that the nitrogen is connected to a carbon (R-N=C=S)
    for match in mol.GetSubstructMatches(isothiocyanate_pattern):
        nitrogen_idx = match[0]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        # Ensure the nitrogen is connected to at least one carbon atom
        if not any(neighbor.GetAtomicNum() == 6 for neighbor in nitrogen_atom.GetNeighbors()):
            return False, "Nitrogen in isothiocyanate group must be connected to a carbon"

    # Additional check to ensure the molecule is not part of a larger, more complex structure
    # For example, exclude molecules with multiple functional groups or complex ring systems
    # that would classify them differently.
    # Here, we check if the molecule has more than one functional group or complex rings.
    # This is a heuristic and may need refinement based on specific cases.
    functional_groups = Chem.MolFromSmarts("[*]!@[N]=[C]=[S]")
    if len(mol.GetSubstructMatches(functional_groups)) > 1:
        return False, "Molecule has multiple functional groups or complex structure"

    return True, "Contains the isothiocyanate group (-N=C=S) and is not part of a larger, more complex structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52221',
                          'name': 'isothiocyanate',
                          'definition': 'An organosulfur compound with the general formula R-N=C=S.',
                          'parents': ['CHEBI:48375', 'CHEBI:48376']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}