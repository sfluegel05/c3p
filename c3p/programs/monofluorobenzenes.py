"""
Classifies: CHEBI:83575 monofluorobenzenes
"""
from rdkit import Chem

def is_monofluorobenzenes(smiles: str):
    """
    Determines if a molecule is a monofluorobenzene (a benzene ring carrying a single fluorine substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monofluorobenzene, False otherwise
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

    # Check if there is exactly one fluorine attached to the benzene ring
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        fluorine_count = 0
        other_atoms = 0
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms:
                    if neighbor.GetSymbol() == 'F':
                        fluorine_count += 1
                    else:
                        other_atoms += 1

        if fluorine_count == 1:
            return True, "Monofluorobenzene found"
        elif fluorine_count > 1:
            return False, "More than one fluorine substituent found"

    return False, "No fluorine substituent found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83575',
                          'name': 'monofluorobenzenes',
                          'definition': 'Any member of the class of '
                                        'fluorobenzenes containing a mono- or '
                                        'poly-substituted benzene ring '
                                        'carrying a single fluorine '
                                        'substitutent.',
                          'parents': ['CHEBI:35496']},
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
    'num_true_positives': 19,
    'num_false_positives': 1,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.95,
    'recall': 1.0,
    'f1': 0.9743589743589743,
    'accuracy': None}