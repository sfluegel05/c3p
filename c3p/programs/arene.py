"""
Classifies: CHEBI:33658 arene
"""
from rdkit import Chem

def is_arene(smiles: str):
    """
    Determines if a molecule is an arene (monocyclic or polycyclic aromatic hydrocarbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one aromatic ring
    aromatic_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetIsAromatic() for atom in atoms):
            aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic rings found"

    return True, "Monocyclic or polycyclic aromatic hydrocarbon (arene)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33658',
                          'name': 'arene',
                          'definition': 'Any monocyclic or polycyclic aromatic '
                                        'hydrocarbon.',
                          'parents': ['CHEBI:33659', 'CHEBI:33663']},
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
    'num_true_positives': 39,
    'num_false_positives': 10,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.7959183673469388,
    'recall': 1.0,
    'f1': 0.8863636363636364,
    'accuracy': None}