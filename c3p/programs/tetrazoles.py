"""
Classifies: CHEBI:35689 tetrazoles
"""
from rdkit import Chem

def is_tetrazoles(smiles: str):
    """
    Determines if a molecule is a tetrazole (a five-membered heterocyclic aromatic skeleton containing four N atoms and one C atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrazole, False otherwise
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

    # Find all 5-membered rings
    five_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]

    # Check for tetrazole characteristics in the 5-membered rings
    for ring in five_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if sum(1 for atom in atoms if atom.GetSymbol() == 'N') == 4 and sum(1 for atom in atoms if atom.GetSymbol() == 'C') == 1:
            if all(atom.GetIsAromatic() for atom in atoms):
                return True, "Tetrazole ring found"

    return False, "No tetrazole ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35689',
                          'name': 'tetrazoles',
                          'definition': 'An azole in which the five-membered '
                                        'heterocyclic aromatic skeleton '
                                        'contains four N atoms and one C atom.',
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
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}