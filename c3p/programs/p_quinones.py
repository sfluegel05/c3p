"""
Classifies: CHEBI:25830 p-quinones
"""
from rdkit import Chem

def is_p_quinones(smiles: str):
    """
    Determines if a molecule is a p-quinone (a quinone with two oxo groups para to each other on a 6-membered ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a p-quinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    if not six_membered_rings:
        return False, "No 6-membered rings found"

    # Check for para-quinone structure in 6-membered rings
    for ring in six_membered_rings:
        oxo_positions = []
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if neighbor.GetSymbol() == 'O' and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                        oxo_positions.append(atom_idx)

        if len(oxo_positions) == 2:
            # Check if the oxo groups are para to each other
            pos1, pos2 = sorted(oxo_positions)
            if (pos2 - pos1) % 6 == 3:
                return True, "Molecule is a p-quinone"

    return False, "No para-quinone structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25830',
                          'name': 'p-quinones',
                          'definition': 'A quinone in which the two oxo groups '
                                        'of the quinone are located para to '
                                        'each other on the 6-membered '
                                        'quinonoid ring.',
                          'parents': ['CHEBI:36141']},
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
    'num_true_positives': 12,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 28,
    'precision': 0.7058823529411765,
    'recall': 0.3,
    'f1': 0.42105263157894735,
    'accuracy': None}