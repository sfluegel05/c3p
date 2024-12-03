"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol (cyclic hemiacetal).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a cyclic structure
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size == 5 or size == 6 for size in ring_sizes):
        return False, "No 5 or 6-membered rings found"

    # Check for the presence of a hemiacetal group
    hemiacetal_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 2:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                carbon1 = neighbors[0]
                carbon2 = neighbors[1]
                if carbon1.GetSymbol() == 'C' and carbon2.GetSymbol() == 'C':
                    carbon1_neighbors = [n.GetSymbol() for n in carbon1.GetNeighbors()]
                    carbon2_neighbors = [n.GetSymbol() for n in carbon2.GetNeighbors()]
                    if 'O' in carbon1_neighbors or 'O' in carbon2_neighbors:
                        hemiacetal_found = True
                        break

    if hemiacetal_found:
        return True, "Cyclic hemiacetal group found"
    else:
        return False, "No cyclic hemiacetal group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38131',
                          'name': 'lactol',
                          'definition': 'Cyclic hemiacetals formed by '
                                        'intramolecular addition of a hydroxy '
                                        'group to an aldehydic or ketonic '
                                        'carbonyl group. They are thus '
                                        '1-oxacycloalkan-2-ols or unsaturated '
                                        'analogues.',
                          'parents': ['CHEBI:5653']},
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
    'num_true_positives': 13,
    'num_false_positives': 6,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.6842105263157895,
    'recall': 1.0,
    'f1': 0.8125000000000001,
    'accuracy': None}