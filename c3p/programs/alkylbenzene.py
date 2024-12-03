"""
Classifies: CHEBI:38976 alkylbenzene
"""
from rdkit import Chem

def is_alkylbenzene(smiles: str):
    """
    Determines if a molecule is an alkylbenzene (benzene substituted with one or more alkyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkylbenzene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered aromatic ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check if all atoms in the aromatic ring are carbon
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(atom.GetSymbol() == 'C' for atom in atoms):
            return False, "Ring contains non-carbon atoms"

    # Check for alkyl substituents
    ring_atoms = set(aromatic_rings[0])
    alkyl_substituents = []

    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in ring_atoms:
                # Check if the substituent is an alkyl group
                if neighbor.GetSymbol() == 'C':
                    if all(n.GetSymbol() in ['C', 'H'] for n in neighbor.GetNeighbors()):
                        alkyl_substituents.append(neighbor.GetSymbol())

    if len(alkyl_substituents) > 0:
        return True, f"Alkylbenzene with alkyl groups: {', '.join(set(alkyl_substituents))}"
    else:
        return False, "No alkyl groups found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38976',
                          'name': 'alkylbenzene',
                          'definition': 'A  monocyclic arene that is benzene '
                                        'substituted with one or more alkyl '
                                        'groups.',
                          'parents': ['CHEBI:22712', 'CHEBI:33847']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 22,
    'num_false_positives': 7,
    'num_true_negatives': 13,
    'num_false_negatives': 6,
    'precision': 0.7586206896551724,
    'recall': 0.7857142857142857,
    'f1': 0.7719298245614034,
    'accuracy': None}