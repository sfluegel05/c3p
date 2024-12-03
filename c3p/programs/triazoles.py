"""
Classifies: CHEBI:35727 triazoles
"""
from rdkit import Chem

def is_triazoles(smiles: str):
    """
    Determines if a molecule is a triazole (an azole with three nitrogen atoms and two carbon atoms in a five-membered ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all 5-membered rings with exactly three nitrogen atoms and two carbon atoms
    triazole_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if sum(1 for atom in atoms if atom.GetSymbol() == 'N') == 3 and sum(1 for atom in atoms if atom.GetSymbol() == 'C') == 2:
                triazole_rings.append(ring)

    if not triazole_rings:
        return False, "No 5-membered rings with three nitrogen atoms and two carbon atoms found"

    return True, "Triazole ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35727',
                          'name': 'triazoles',
                          'definition': 'An azole in which the five-membered '
                                        'heterocyclic aromatic skeleton '
                                        'contains three N atoms and two C '
                                        'atoms.',
                          'parents': ['CHEBI:68452']},
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
    'num_true_positives': 43,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}