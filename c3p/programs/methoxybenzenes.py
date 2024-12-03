"""
Classifies: CHEBI:51683 methoxybenzenes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_methoxybenzenes(smiles: str):
    """
    Determines if a molecule is a methoxybenzene (aromatic ether with benzene substituted with methoxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methoxybenzene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered aromatic ring
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

    # Check if the aromatic rings have methoxy groups attached
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        methoxy_groups = []

        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 0:
                    for second_neighbor in neighbor.GetNeighbors():
                        if second_neighbor.GetSymbol() == 'C' and second_neighbor.GetTotalNumHs() == 3:
                            methoxy_groups.append((atom_idx, neighbor.GetIdx(), second_neighbor.GetIdx()))

        if methoxy_groups:
            return True, f"Methoxybenzene with methoxy groups at positions: {', '.join([str(m[0]) for m in methoxy_groups])}"

    return False, "No methoxy groups found on the aromatic ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51683',
                          'name': 'methoxybenzenes',
                          'definition': 'Any aromatic ether that consists of a '
                                        'benzene skeleton substituted with one '
                                        'or more methoxy groups.',
                          'parents': ['CHEBI:22712', 'CHEBI:35618']},
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
    'num_true_positives': 170,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}