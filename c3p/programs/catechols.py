"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol (contains an o-diphenol component).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered aromatic ring
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check for o-diphenol (ortho-dihydroxy) groups in the aromatic rings
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        for i, atom in enumerate(atoms):
            if atom.GetSymbol() == 'C' and atom.GetIsAromatic():
                neighbors = atom.GetNeighbors()
                hydroxyl_neighbors = [n for n in neighbors if n.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1]
                if len(hydroxyl_neighbors) == 1:
                    for j in range(i+1, len(atoms)):
                        next_atom = atoms[j]
                        if next_atom.GetSymbol() == 'C' and next_atom.GetIsAromatic():
                            next_neighbors = next_atom.GetNeighbors()
                            next_hydroxyl_neighbors = [n for n in next_neighbors if n.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(next_atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1]
                            if len(next_hydroxyl_neighbors) == 1:
                                if abs(i - j) == 1 or abs(i - j) == 5:
                                    return True, "Contains o-diphenol (catechol) component"

    return False, "No o-diphenol (catechol) component found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33566',
                          'name': 'catechols',
                          'definition': 'Any compound containing an o-diphenol '
                                        'component.',
                          'parents': ['CHEBI:33570']},
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
    'num_true_positives': 39,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9512195121951219,
    'recall': 1.0,
    'f1': 0.975,
    'accuracy': None}