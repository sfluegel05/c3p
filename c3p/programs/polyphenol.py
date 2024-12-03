"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Find all aromatic rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetIsAromatic() for atom in atoms):
            aromatic_rings.append(ring)

    if len(aromatic_rings) < 2:
        return False, "Less than 2 aromatic rings found"

    # Check if each aromatic ring has at least one hydroxyl group
    ring_hydroxy_count = 0
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        has_hydroxy = False
        for atom in atoms:
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() > 0:
                        has_hydroxy = True
                        break
            if has_hydroxy:
                break
        if has_hydroxy:
            ring_hydroxy_count += 1

    if ring_hydroxy_count < 2:
        return False, "Not all aromatic rings have a hydroxyl group"

    return True, "Molecule is a polyphenol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26195',
                          'name': 'polyphenol',
                          'definition': 'Members of the class of phenols that '
                                        'contain 2 or more benzene rings each '
                                        'of which is substituted by at least '
                                        'one hydroxy group.',
                          'parents': ['CHEBI:33853']},
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
    'num_true_positives': 73,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 12,
    'precision': 0.9605263157894737,
    'recall': 0.8588235294117647,
    'f1': 0.906832298136646,
    'accuracy': None}