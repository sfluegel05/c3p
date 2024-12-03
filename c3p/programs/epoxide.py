"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide (cyclic ether with a 3-membered ring containing an oxygen atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
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

    # Find all 3-membered rings containing an oxygen atom
    epoxide_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 3:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                epoxide_rings.append(ring)

    if not epoxide_rings:
        return False, "No 3-membered rings containing an oxygen atom found"

    return True, "Epoxide detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32955',
                          'name': 'epoxide',
                          'definition': 'Any cyclic ether in which the oxygen '
                                        'atom forms part of a 3-membered ring.',
                          'parents': ['CHEBI:37407']},
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
    'num_true_positives': 82,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9879518072289156,
    'f1': 0.993939393939394,
    'accuracy': None}