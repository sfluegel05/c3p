"""
Classifies: CHEBI:78802 diarylheptanoid
"""
from rdkit import Chem

def is_diarylheptanoid(smiles: str):
    """
    Determines if a molecule is a diarylheptanoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diarylheptanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least two phenyl rings
    phenyl_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() and atom.GetSymbol() == 'C' for atom in atoms):
                phenyl_rings.append(ring)

    if len(phenyl_rings) < 2:
        return False, "Less than two phenyl rings found"

    # Check for heptane chain connecting the two phenyl rings
    def find_path_length(atom1, atom2):
        queue = [(atom1, 0)]
        visited = set()
        while queue:
            current_atom, length = queue.pop(0)
            if current_atom == atom2:
                return length
            visited.add(current_atom)
            for neighbor in mol.GetAtomWithIdx(current_atom).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetSymbol() == 'C':
                    queue.append((neighbor_idx, length + 1))
        return float('inf')

    for i in range(len(phenyl_rings)):
        for j in range(i + 1, len(phenyl_rings)):
            for atom1 in phenyl_rings[i]:
                for atom2 in phenyl_rings[j]:
                    if 7 <= find_path_length(atom1, atom2) <= 10:  # Adjusted to include substituted heptane chains
                        return True, "Molecule is a diarylheptanoid"

    return False, "No heptane chain connecting the two phenyl rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:78802',
                          'name': 'diarylheptanoid',
                          'definition': 'A family of plant metabolites with a '
                                        'common 1,7-diphenylheptane structural '
                                        'skeleton, carrying various '
                                        'substituents. They are mainly '
                                        'distributed in the roots, rhizomes '
                                        'and bark of Alpinia, Zingiber, '
                                        'Curcuma and Alnus species.',
                          'parents': ['CHEBI:33822', 'CHEBI:33836']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 21,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9130434782608695,
    'recall': 1.0,
    'f1': 0.9545454545454545,
    'accuracy': None}