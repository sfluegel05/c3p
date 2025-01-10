"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:52214 spiroketal
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible spiroketal pattern: a carbon atom bonded to two oxygen atoms, each in a separate ring
    spiroketal_pattern = Chem.MolFromSmarts("[CX4]([OX2][*;r])[OX2][*;r]")
    
    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    
    if not matches:
        return False, "No spiroketal pattern found"
    
    # Check if the ketal carbon is the only common atom between the two rings
    for match in matches:
        ketal_carbon_idx = match[0]
        oxygen1_idx = match[1]
        oxygen2_idx = match[2]
        
        # Get the rings that contain the oxygens
        rings = mol.GetRingInfo().AtomRings()
        ring1 = None
        ring2 = None
        
        for ring in rings:
            if oxygen1_idx in ring:
                ring1 = ring
            if oxygen2_idx in ring:
                ring2 = ring
        
        if ring1 is None or ring2 is None:
            continue
        
        # Check if the ketal carbon is the only common atom between the two rings
        common_atoms = set(ring1).intersection(set(ring2))
        if len(common_atoms) == 1 and ketal_carbon_idx in common_atoms:
            return True, "Contains a spiroketal structure with a ketal carbon as the only common atom between two rings"
    
    return False, "No spiroketal structure found with a ketal carbon as the only common atom between two rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52214',
                          'name': 'spiroketal',
                          'definition': 'A cyclic ketal in which the ketal '
                                        'carbon is the only common atom of '
                                        'two rings.',
                          'parents': ['CHEBI:52213', 'CHEBI:52215']},
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