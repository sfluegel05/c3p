"""
Classifies: CHEBI:71061 octasaccharide derivative
"""
from rdkit import Chem

def is_octasaccharide_derivative(smiles: str):
    """
    Determines if a molecule is an octasaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octasaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify sugar units by checking for typical sugar rings (5 or 6 membered rings with oxygen)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_units = 0

    for ring in rings:
        if len(ring) == 5 or len(ring) == 6:
            if any(mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'O' for atom_idx in ring):
                sugar_units += 1

    if sugar_units < 8:
        return False, f"Only {sugar_units} sugar units found, less than 8"

    # Check for derivative groups
    derivative_groups = ['N', 'S', 'P']
    has_derivative = any(atom.GetSymbol() in derivative_groups for atom in mol.GetAtoms())

    if not has_derivative:
        return False, "No derivative groups found"

    return True, "Molecule is an octasaccharide derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:71061',
                          'name': 'octasaccharide derivative',
                          'definition': 'An oligosaccharide derivative that is '
                                        'formally obtained from an '
                                        'octasaccharide.',
                          'parents': ['CHEBI:63563']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 1,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.9090909090909091,
    'recall': 1.0,
    'f1': 0.9523809523809523,
    'accuracy': None}