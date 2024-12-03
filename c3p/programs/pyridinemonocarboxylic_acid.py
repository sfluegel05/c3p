"""
Classifies: CHEBI:26420 pyridinemonocarboxylic acid
"""
from rdkit import Chem

def is_pyridinemonocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a pyridinemonocarboxylic acid (a monocarboxylic acid in which the carboxy group is attached to a pyridine or substituted pyridine ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyridinemonocarboxylic acid, False otherwise
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

    # Find all pyridine rings (aromatic 6-membered rings with one nitrogen)
    pyridine_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if sum(1 for atom in atoms if atom.GetSymbol() == 'N') == 1 and all(atom.GetIsAromatic() for atom in atoms):
                pyridine_rings.append(ring)

    if not pyridine_rings:
        return False, "No pyridine rings found"

    # Check for carboxylic acid group attached to the pyridine ring
    carboxylic_acid_found = False
    for ring in pyridine_rings:
        ring_atoms = set(ring)
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    # Check if the carbon has two oxygen neighbors (one double-bonded and one single-bonded)
                    oxygen_neighbors = [n for n in neighbor.GetNeighbors() if n.GetSymbol() == 'O']
                    if len(oxygen_neighbors) == 2:
                        if any(bond.GetBondTypeAsDouble() == 2 for bond in neighbor.GetBonds() if bond.GetOtherAtom(neighbor).GetSymbol() == 'O'):
                            carboxylic_acid_found = True
                            break
            if carboxylic_acid_found:
                break
        if carboxylic_acid_found:
            break

    if carboxylic_acid_found:
        return True, "Pyridinemonocarboxylic acid found"
    else:
        return False, "No carboxylic acid group attached to pyridine ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26420',
                          'name': 'pyridinemonocarboxylic acid',
                          'definition': 'A monocarboxylic acid in which the '
                                        'carboxy group is attached to a '
                                        'pyridine (or substituted pyridine) '
                                        'ring.',
                          'parents': [   'CHEBI:25384',
                                         'CHEBI:26421',
                                         'CHEBI:33859']},
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
    'num_true_positives': 7,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.7,
    'f1': 0.8235294117647058,
    'accuracy': None}