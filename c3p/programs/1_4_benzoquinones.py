"""
Classifies: CHEBI:132124 1,4-benzoquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_1_4_benzoquinones(smiles: str):
    """
    Determines if a molecule is a 1,4-benzoquinone or its C-substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,4-benzoquinone or its C-substituted derivatives, False otherwise
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

    # Find all 1,4-benzoquinone-like rings
    benzoquinone_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            carbonyl_positions = []
            for atom in atoms:
                if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            carbonyl_positions.append(atom.GetIdx())
            if len(carbonyl_positions) == 2:
                # Ensure the carbonyl groups are at positions 1 and 4
                pos_diff = abs(ring.index(carbonyl_positions[0]) - ring.index(carbonyl_positions[1]))
                if pos_diff == 3:
                    benzoquinone_rings.append(ring)

    if not benzoquinone_rings:
        return False, "No 1,4-benzoquinone-like rings found"

    return True, "1,4-benzoquinone or its C-substituted derivatives found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:132124',
                          'name': '1,4-benzoquinones',
                          'definition': 'Any member of the class of '
                                        'benzoquinones that is '
                                        '1,4-benzoquinone or its C-substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:22729', 'CHEBI:25830']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 7,
    'num_true_negatives': 3,
    'num_false_negatives': 0,
    'precision': 0.5882352941176471,
    'recall': 1.0,
    'f1': 0.7407407407407407,
    'accuracy': None}