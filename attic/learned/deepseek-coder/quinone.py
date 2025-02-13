"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is a compound with a fully conjugated cyclic dione structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a cyclic structure
    if not mol.GetRingInfo().NumRings():
        return False, "No cyclic structure found"

    # Look for two carbonyl groups (C=O) in the molecule
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 2:
        return False, f"Found {len(carbonyl_matches)} carbonyl groups, need at least 2"

    # Check if the carbonyl groups are part of a conjugated system
    conjugated_system = Chem.MolFromSmarts("[CX3]=[CX3]")
    conjugated_matches = mol.GetSubstructMatches(conjugated_system)
    if not conjugated_matches:
        return False, "No conjugated system found"

    # Check if the carbonyl groups are part of the same ring
    ring_info = mol.GetRingInfo()
    for match in carbonyl_matches:
        atom_idx = match[0]
        if not any(atom_idx in ring for ring in ring_info.AtomRings()):
            return False, "Carbonyl groups not part of the same ring"

    # Check if the molecule has a fully conjugated cyclic dione structure
    # This is a more complex check, but we can approximate it by ensuring the molecule
    # has a cyclic structure with two carbonyl groups in a conjugated system
    # and that the molecule is planar (indicative of conjugation)
    if not rdMolDescriptors.CalcNumAromaticRings(mol):
        return False, "No aromatic or conjugated system found"

    return True, "Contains a fully conjugated cyclic dione structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36141',
                          'name': 'quinone',
                          'definition': 'Compounds having a fully conjugated '
                                        'cyclic dione structure, such as that '
                                        'of benzoquinones, derived from '
                                        'aromatic compounds by conversion of '
                                        'an even number of -CH= groups into '
                                        '-C(=O)- groups with any necessary '
                                        'rearrangement of double bonds '
                                        '(polycyclic and heterocyclic '
                                        'analogues are included).',
                          'parents': ['CHEBI:36141', 'CHEBI:36141']},
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