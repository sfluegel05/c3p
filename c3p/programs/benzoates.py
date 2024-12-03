"""
Classifies: CHEBI:22718 benzoates
"""
from rdkit import Chem

def is_benzoates(smiles: str):
    """
    Determines if a molecule is a benzoate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzoate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
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

    # Check if the molecule contains a carboxylate group ([O-]C=O)
    carboxylate_group = Chem.MolFromSmarts('[O-]C=O')
    if not mol.HasSubstructMatch(carboxylate_group):
        return False, "No carboxylate group found"

    # Check if the carboxylate group is attached to the aromatic ring
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms:
                    if mol.GetSubstructMatch(carboxylate_group, useChirality=True):
                        return True, "Molecule is a benzoate"

    return False, "Carboxylate group is not attached to the aromatic ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22718',
                          'name': 'benzoates',
                          'definition': 'A monocarboxylic acid anion obtained '
                                        'by deprotonation of the carboxy group '
                                        'of any benzoic acid.',
                          'parents': ['CHEBI:35757', 'CHEBI:91007']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 14,
    'num_false_positives': 4,
    'num_true_negatives': 12,
    'num_false_negatives': 2,
    'precision': 0.7777777777777778,
    'recall': 0.875,
    'f1': 0.823529411764706,
    'accuracy': None}