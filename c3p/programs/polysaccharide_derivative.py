"""
Classifies: CHEBI:65212 polysaccharide derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide_derivative(smiles: str):
    """
    Determines if a molecule is a polysaccharide derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify polysaccharide characteristics
    # Polysaccharides are long chains of monosaccharide units bound together by glycosidic linkages
    # We will look for repeating saccharide units and glycosidic bonds (O-C-O)

    # Check for the presence of multiple saccharide units
    saccharide_count = 0
    for ring in mol.GetRingInfo().AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if len(atoms) == 5 or len(atoms) == 6:  # Typical ring sizes in saccharides
            if all(atom.GetSymbol() == 'C' or atom.GetSymbol() == 'O' for atom in atoms):
                saccharide_count += 1

    if saccharide_count < 2:
        return False, "Less than two saccharide units found"

    # Check for glycosidic linkages (O-C-O)
    glycosidic_linkages = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2 and all(neigh.GetSymbol() == 'C' for neigh in neighbors):
                glycosidic_linkages += 1

    if glycosidic_linkages < 1:
        return False, "No glycosidic linkages (O-C-O) found"

    return True, "Polysaccharide derivative identified"

# Example usage:
# smiles = "O1C(OC(CC(O)=O)C)C(O)C(O)C(O)C1C(O)=O"
# print(is_polysaccharide_derivative(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65212',
                          'name': 'polysaccharide derivative',
                          'definition': 'A carbohydrate derivative that is any '
                                        'derivative of a polysaccharide.',
                          'parents': [   'CHEBI:167559',
                                         'CHEBI:33694',
                                         'CHEBI:63299']},
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
    'num_true_positives': 28,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 18,
    'precision': 0.5833333333333334,
    'recall': 0.6086956521739131,
    'f1': 0.5957446808510638,
    'accuracy': None}