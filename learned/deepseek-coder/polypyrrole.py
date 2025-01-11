"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: CHEBI:38834 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define multiple patterns to catch different pyrrole configurations
    patterns = [
        # Standard pyrrole ring
        Chem.MolFromSmarts("[n;H0,H1]1[c;H1,H0][c;H1,H0][c;H1,H0][c;H1,H0]1"),
        # Pyrrole with double bonds
        Chem.MolFromSmarts("[n;H0,H1]1[c;H1,H0]=[c;H1,H0][c;H1,H0][c;H1,H0]1"),
        # Pyrrole with substitutions
        Chem.MolFromSmarts("[n;H0,H1]1[c;H1,H0][c;H1,H0][c;H1,H0][c;H1,H0]1(-*)"),
        # Fused pyrrole systems
        Chem.MolFromSmarts("[n;H0,H1]1[c;H1,H0][c;H1,H0][c;H1,H0][c;H1,H0]1~*"),
        # Conjugated pyrrole-like structures
        Chem.MolFromSmarts("[n;H0,H1]1[c;H1,H0][c;H1,H0][c;H1,H0][c;H1,H0]1~[c,n]"),
    ]
    
    # Count all unique pyrrole units using all patterns
    pyrrole_units = set()
    for pattern in patterns:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # Store the atom indices to avoid double counting
            pyrrole_units.add(tuple(sorted(match)))
    
    # Check if we have at least two distinct pyrrole units
    if len(pyrrole_units) >= 2:
        return True, f"Contains {len(pyrrole_units)} pyrrole units"
    
    # If not, try to find pyrrole-like structures in conjugated systems
    conjugated_pattern = Chem.MolFromSmarts("[n;H0,H1]~[c,n]~[c,n]~[c,n]~[c,n]")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) >= 2:
        return True, "Contains conjugated pyrrole-like structures"
    
    return False, f"Found {len(pyrrole_units)} pyrrole units, need at least 2"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38834',
                          'name': 'polypyrrole',
                          'definition': 'A compound composed of two or more pyrrole units.',
                          'parents': ['CHEBI:38834']},
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