"""
Classifies: CHEBI:20857 C-glycosyl compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_C_glycosyl_compound(smiles: str):
    """
    Determines if a molecule is a C-glycosyl compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a C-glycosyl compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify all glycosidic rings (oxane rings)
    rings = mol.GetRingInfo()
    oxane_rings = [ring for ring in rings.AtomRings() if len(ring) == 6 and any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in ring)]

    if not oxane_rings:
        return False, "No oxane rings found"

    # Check for C-glycosidic bond (C-C bond between glycosidic ring and another carbon)
    for ring in oxane_rings:
        ring_atoms = set(ring)
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in ring_atoms and neighbor.GetSymbol() == 'C':
                        # Ensure the neighbor carbon is not part of a glycosidic hydroxy group
                        if not any(neigh.GetSymbol() == 'O' for neigh in neighbor.GetNeighbors()):
                            return True, "C-glycosidic bond (C-C bond) found"

    return False, "No C-glycosidic bond (C-C bond) found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:20857',
                          'name': 'C-glycosyl compound',
                          'definition': 'A glycosyl compound arising formally '
                                        'from the elimination of water from a '
                                        'glycosidic hydroxy group and an H '
                                        'atom bound to a carbon atom, thus '
                                        'creating a C-C bond.',
                          'parents': ['CHEBI:63161', 'CHEBI:63299']},
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
    'num_true_positives': 214,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 10,
    'precision': 0.9727272727272728,
    'recall': 0.9553571428571429,
    'f1': 0.963963963963964,
    'accuracy': None}