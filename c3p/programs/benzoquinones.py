"""
Classifies: CHEBI:22729 benzoquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_benzoquinones(smiles: str):
    """
    Determines if a molecule is a benzoquinone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzoquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all 6-membered rings
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]

    # Check for 1,4-benzoquinone structure
    for ring in six_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        carbonyl_groups = 0
        for i, atom in enumerate(atoms):
            if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        carbonyl_groups += 1
                        break

        if carbonyl_groups == 2:
            # Check for catechol or hydroquinone oxidation pattern
            if any(atom.GetSymbol() == 'C' and atom.GetDegree() == 3 and 
                   any(neigh.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neigh.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE 
                       for neigh in atom.GetNeighbors()) for atom in atoms):
                return True, "1,4-benzoquinone structure found"

    return False, "No 1,4-benzoquinone structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22729',
                          'name': 'benzoquinones',
                          'definition': 'Any quinone resulting from the formal '
                                        'oxidation of catechol, hydroquinone, '
                                        'or their C-substituted derivatives.',
                          'parents': ['CHEBI:36141']},
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
    'num_true_positives': 14,
    'num_false_positives': 5,
    'num_true_negatives': 12,
    'num_false_negatives': 3,
    'precision': 0.7368421052631579,
    'recall': 0.8235294117647058,
    'f1': 0.7777777777777778,
    'accuracy': None}