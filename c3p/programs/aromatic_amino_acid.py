"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid (an amino acid whose structure includes an aromatic ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                carboxylic_acid = True
                break
    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for presence of amino group (NH2 or NH)
    amino_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetDegree() in [1, 2]:
            amino_group = True
            break
    if not amino_group:
        return False, "No amino group found"

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

    return True, "Aromatic amino acid found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33856',
                          'name': 'aromatic amino acid',
                          'definition': 'An amino acid whose structure '
                                        'includes an aromatic ring.',
                          'parents': ['CHEBI:33659', 'CHEBI:33709']},
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
    'num_true_positives': 11,
    'num_false_positives': 3,
    'num_true_negatives': 8,
    'num_false_negatives': 0,
    'precision': 0.7857142857142857,
    'recall': 1.0,
    'f1': 0.88,
    'accuracy': None}