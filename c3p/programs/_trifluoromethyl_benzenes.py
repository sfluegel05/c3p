"""
Classifies: CHEBI:83565 (trifluoromethyl)benzenes
"""
from rdkit import Chem

def is_trifluoromethyl_benzenes(smiles: str):
    """
    Determines if a molecule is a (trifluoromethyl)benzene or its derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a (trifluoromethyl)benzene or its derivative, False otherwise
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

    # Check if any aromatic ring has a trifluoromethyl group (CF3)
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 0:
                    neighbor_neighbors = neighbor.GetNeighbors()
                    if len(neighbor_neighbors) == 3 and all(nn.GetSymbol() == 'F' for nn in neighbor_neighbors):
                        return True, "Contains a trifluoromethyl group attached to an aromatic ring"

    return False, "No trifluoromethyl group attached to an aromatic ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83565',
                          'name': '(trifluoromethyl)benzenes',
                          'definition': 'An organofluorine compound that is '
                                        '(trifluoromethyl)benzene and '
                                        'derivatives arising from substitution '
                                        'of one or more of the phenyl '
                                        'hydrogens.',
                          'parents': ['CHEBI:22712', 'CHEBI:37143']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "name 'is__trifluoromethyl_benzenes' is not defined",
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}