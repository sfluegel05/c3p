"""
Classifies: CHEBI:38932 pyridopyrimidine
"""
from rdkit import Chem

def is_pyridopyrimidine(smiles: str):
    """
    Determines if a molecule is a pyridopyrimidine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyridopyrimidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all 6-membered rings
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]

    # Check for the presence of both pyridine and pyrimidine rings
    pyridine_rings = []
    pyrimidine_rings = []

    for ring in six_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        atom_symbols = [atom.GetSymbol() for atom in atoms]
        num_nitrogens = atom_symbols.count('N')

        if num_nitrogens == 1 and 'N' in atom_symbols:
            pyridine_rings.append(ring)
        elif num_nitrogens == 2 and 'N' in atom_symbols:
            pyrimidine_rings.append(ring)

    if not pyridine_rings:
        return False, "No pyridine rings found"
    if not pyrimidine_rings:
        return False, "No pyrimidine rings found"

    # Check if pyridine and pyrimidine rings are ortho-fused
    for pyridine_ring in pyridine_rings:
        for pyrimidine_ring in pyrimidine_rings:
            shared_atoms = set(pyridine_ring).intersection(set(pyrimidine_ring))
            if len(shared_atoms) >= 2:
                return True, "Pyridopyrimidine detected"

    return False, "Pyridine and pyrimidine rings are not ortho-fused"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38932',
                          'name': 'pyridopyrimidine',
                          'definition': 'Any organic heterobicyclic compound '
                                        'consisting of a pyridine ring '
                                        'ortho-fused at any position to a '
                                        'pyrimidine ring.',
                          'parents': ['CHEBI:27171']},
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
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}