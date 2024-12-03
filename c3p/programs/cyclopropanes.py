"""
Classifies: CHEBI:51454 cyclopropanes
"""
from rdkit import Chem

def is_cyclopropanes(smiles: str):
    """
    Determines if a molecule is a cyclopropane or a derivative formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclopropane or derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 3-membered ring
    if not any(len(ring) == 3 for ring in rings.AtomRings()):
        return False, "No 3-membered rings found"

    # Find all 3-membered rings
    cyclopropane_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 3:
            cyclopropane_rings.append(ring)

    if not cyclopropane_rings:
        return False, "No cyclopropane rings found"

    # Check if the 3-membered ring is composed of carbon atoms
    for ring in cyclopropane_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(atom.GetSymbol() == 'C' for atom in atoms):
            return False, "3-membered ring contains non-carbon atoms"

    # Check substituents
    ring_atoms = set(cyclopropane_rings[0])
    substituents = []

    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in ring_atoms:
                substituents.append(neighbor.GetSymbol())

    if len(substituents) > 0:
        return True, f"Cyclopropane derivative with substituents: {', '.join(set(substituents))}"
    else:
        return True, "Unsubstituted cyclopropane"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51454',
                          'name': 'cyclopropanes',
                          'definition': 'Cyclopropane and its derivatives '
                                        'formed by substitution.',
                          'parents': ['CHEBI:33598']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[01:08:54] SMILES Parse Error: syntax error while parsing: '
             'CO[C@@H]1CC(=O)Nc2cc(O)cc(CC\\C=C(C)/[C@H](O)[C@@H](C)[C@H](C\\C=C\\C=C\\C=C\x01)OC(=O)C1(CC1)NC(=O)C1=CCCCC1)c2O\n'
             '[01:08:54] SMILES Parse Error: Failed parsing SMILES '
             "'CO[C@@H]1CC(=O)Nc2cc(O)cc(CC\\C=C(C)/[C@H](O)[C@@H](C)[C@H](C\\C=C\\C=C\\C=C\x01)OC(=O)C1(CC1)NC(=O)C1=CCCCC1)c2O' "
             'for input: '
             "'CO[C@@H]1CC(=O)Nc2cc(O)cc(CC\\C=C(C)/[C@H](O)[C@@H](C)[C@H](C\\C=C\\C=C\\C=C\x01)OC(=O)C1(CC1)NC(=O)C1=CCCCC1)c2O'\n",
    'stdout': '',
    'num_true_positives': 20,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9523809523809523,
    'f1': 0.975609756097561,
    'accuracy': None}