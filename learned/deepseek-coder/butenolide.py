"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: CHEBI:37741 butenolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone with a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core butenolide pattern (2-furanone skeleton)
    butenolide_pattern = Chem.MolFromSmarts("[O;R]=[C;R][C;R]=[C;R][C;R]=O")
    
    # Check for the core pattern
    if not mol.HasSubstructMatch(butenolide_pattern):
        return False, "No 2-furanone skeleton found"

    # Verify it's a 5-membered ring
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not any(len(ring) == 5 for ring in atom_rings):
        return False, "No 5-membered ring found"

    # Count ester groups (should have at least one)
    ester_pattern = Chem.MolFromSmarts("[O;R]=[C;R][O;R]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester group found"

    # Check for double bond in the ring
    double_bond_pattern = Chem.MolFromSmarts("[C;R]=[C;R]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 1:
        return False, "No double bond found in the ring"

    return True, "Contains 2-furanone skeleton with gamma-lactone structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:37741',
        'name': 'butenolide',
        'definition': 'A gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.',
        'parents': ['CHEBI:37671', 'CHEBI:51086']
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}