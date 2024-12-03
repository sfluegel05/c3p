"""
Classifies: CHEBI:35990 bridged compound
"""
from rdkit import Chem

def is_bridged_compound(smiles: str):
    """
    Determines if a molecule is a bridged compound (contains more than one ring with at least two common atoms that are not adjacent to each other).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bridged compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()
    atom_rings = rings.AtomRings()

    if len(atom_rings) < 2:
        return False, "Less than two rings found"

    # Find common atoms between rings
    common_atoms = {}
    for i, ring1 in enumerate(atom_rings):
        for j, ring2 in enumerate(atom_rings):
            if i >= j:
                continue
            common = set(ring1).intersection(set(ring2))
            if len(common) >= 2:
                for atom in common:
                    if atom not in common_atoms:
                        common_atoms[atom] = []
                    common_atoms[atom].append((i, j))

    if not common_atoms:
        return False, "No common atoms found between rings"

    # Check if common atoms are bridgehead atoms (not adjacent)
    for atom, ring_pairs in common_atoms.items():
        for ring1, ring2 in ring_pairs:
            if any(abs(ring1.index(atom) - ring1.index(other_atom)) > 1 for other_atom in atom_rings[ring1]):
                if any(abs(ring2.index(atom) - ring2.index(other_atom)) > 1 for other_atom in atom_rings[ring2]):
                    return True, "Bridged compound with bridgehead atoms"

    return False, "No bridgehead atoms found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35990',
                          'name': 'bridged compound',
                          'definition': 'A polycyclic compound that contains '
                                        'more than one ring with at least two '
                                        'common atoms (also known as '
                                        'bridgehead carbons) that are not '
                                        'adjacent to each other.',
                          'parents': ['CHEBI:33635']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "'int' object has no attribute 'index'",
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}