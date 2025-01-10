"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:47778 2,5-diketopiperazine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the piperazine-2,5-dione skeleton pattern
    # The pattern matches a six-membered ring with two nitrogen atoms at positions 1 and 4,
    # and two carbonyl groups at positions 2 and 5.
    piperazine_dione_pattern = Chem.MolFromSmarts("O=C1NCC(=O)NC1")
    
    # Check if the molecule contains the piperazine-2,5-dione skeleton
    if not mol.HasSubstructMatch(piperazine_dione_pattern):
        return False, "No piperazine-2,5-dione skeleton found"

    # Additional checks to ensure the molecule is a 2,5-diketopiperazine
    # Count the number of nitrogen and oxygen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # A 2,5-diketopiperazine must have at least 2 nitrogen and 2 oxygen atoms
    if n_count < 2:
        return False, "Not enough nitrogen atoms (need at least 2)"
    if o_count < 2:
        return False, "Not enough oxygen atoms (need at least 2)"

    # Check for the presence of the piperazine ring
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    piperazine_ring_found = False
    
    for ring in rings:
        if len(ring) == 6:  # Piperazine is a 6-membered ring
            # Check if the ring contains two nitrogen atoms
            nitrogen_count = sum(1 for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7)
            if nitrogen_count == 2:
                piperazine_ring_found = True
                break
    
    if not piperazine_ring_found:
        return False, "No piperazine ring found"

    return True, "Contains a piperazine-2,5-dione skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47778',
                          'name': '2,5-diketopiperazine',
                          'definition': 'Any piperazinone that has a '
                                        'piperazine-2,5-dione skeleton.',
                          'parents': ['CHEBI:47778', 'CHEBI:76579']},
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