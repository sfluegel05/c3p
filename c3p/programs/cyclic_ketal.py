"""
Classifies: CHEBI:59779 cyclic ketal
"""
from rdkit import Chem

def is_cyclic_ketal(smiles: str):
    """
    Determines if a molecule is a cyclic ketal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic ketal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms to find possible ketal carbons
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Check for carbon atom
            neighbors = atom.GetNeighbors()
            oxygen_count = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 8)
            if oxygen_count == 2:
                # Check if both oxygens are part of rings
                ring_info = mol.GetRingInfo()
                atom_idx = atom.GetIdx()
                ring_atoms = set()
                for ring in ring_info.AtomRings():
                    if atom_idx in ring:
                        ring_atoms.update(ring)
                if all(neighbor.GetIdx() in ring_atoms for neighbor in neighbors if neighbor.GetAtomicNum() == 8):
                    return True, "Cyclic ketal structure found"
    
    return False, "No cyclic ketal structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59779',
                          'name': 'cyclic ketal',
                          'definition': 'A ketal in the molecule of which the '
                                        'ketal carbon and one or both oxygen '
                                        'atoms thereon are members of a ring.',
                          'parents': ['CHEBI:38104', 'CHEBI:59777']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 16,
    'num_false_positives': 13,
    'num_true_negatives': 4,
    'num_false_negatives': 1,
    'precision': 0.5517241379310345,
    'recall': 0.9411764705882353,
    'f1': 0.6956521739130435,
    'accuracy': None}