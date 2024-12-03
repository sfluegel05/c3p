"""
Classifies: CHEBI:25481 naphthoquinone
"""
from rdkit import Chem

def is_naphthoquinone(smiles: str):
    """
    Determines if a molecule is a naphthoquinone (a polycyclic aromatic ketone metabolite of naphthalene).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a naphthoquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least two 6-membered rings
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    if len(six_membered_rings) < 2:
        return False, "Less than two 6-membered rings found"

    # Check if the molecule contains ketone groups (C=O)
    ketone_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3 and any(nb.GetSymbol() == 'O' and nb.GetTotalDegree() == 1 for nb in atom.GetNeighbors())]
    if len(ketone_groups) < 2:
        return False, "Less than two ketone groups found"

    # Check for naphthalene core structure
    naphthalene_core = False
    for ring1 in six_membered_rings:
        for ring2 in six_membered_rings:
            if ring1 != ring2 and any(atom in ring1 for atom in ring2):
                naphthalene_core = True
                break
        if naphthalene_core:
            break

    if not naphthalene_core:
        return False, "No naphthalene core structure found"

    return True, "Naphthoquinone structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25481',
                          'name': 'naphthoquinone',
                          'definition': 'A polycyclic aromatic ketone '
                                        'metabolite of naphthalene.',
                          'parents': ['CHEBI:36141']},
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
    'num_true_positives': 25,
    'num_false_positives': 18,
    'num_true_negatives': 2,
    'num_false_negatives': 0,
    'precision': 0.5813953488372093,
    'recall': 1.0,
    'f1': 0.7352941176470588,
    'accuracy': None}