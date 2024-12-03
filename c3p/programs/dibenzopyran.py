"""
Classifies: CHEBI:39203 dibenzopyran
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dibenzopyran(smiles: str):
    """
    Determines if a molecule is a dibenzopyran (a pyran ring fused with two benzene rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dibenzopyran, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Find all 6-membered rings
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    if len(six_membered_rings) < 3:
        return False, "Less than three 6-membered rings found"

    # Check for at least one 6-membered ring with an oxygen atom (pyran ring)
    pyran_rings = [ring for ring in six_membered_rings if any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in ring)]
    if not pyran_rings:
        return False, "No pyran rings found"

    # Check for fusion with two benzene rings
    for pyran_ring in pyran_rings:
        adjacent_rings = []
        for atom_idx in pyran_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                for ring in six_membered_rings:
                    if neighbor_idx in ring and ring not in adjacent_rings and ring != pyran_ring:
                        adjacent_rings.append(ring)
        if len(adjacent_rings) >= 2:
            return True, "Dibenzopyran structure found"

    return False, "No dibenzopyran structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:39203',
                          'name': 'dibenzopyran',
                          'definition': 'Any organic heteropolycyclic compound '
                                        'based on a skeleton consisting of a '
                                        'pyran ring fused with two benzene '
                                        'rings.',
                          'parents': ['CHEBI:38104', 'CHEBI:38166']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 82,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 1,
    'precision': 0.9879518072289156,
    'recall': 0.9879518072289156,
    'f1': 0.9879518072289156,
    'accuracy': None}