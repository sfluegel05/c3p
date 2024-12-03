"""
Classifies: CHEBI:33245 organic fundamental parent
"""
from rdkit import Chem

def is_organic_fundamental_parent(smiles: str):
    """
    Determines if a molecule is an organic fundamental parent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic fundamental parent, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carbon atoms
    has_carbon = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            has_carbon = True
            break

    if not has_carbon:
        return False, "No carbon atoms found"

    # Check for single atom of an element, unbranched chain, monocyclic or polycyclic ring system, or ring/chain system
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()

    if num_atoms == 1:
        return True, "Single atom of an element"
    elif num_rings == 0 and num_bonds + 1 == num_atoms:
        return True, "Unbranched chain"
    elif num_rings > 0:
        return True, "Monocyclic or polycyclic ring system"
    else:
        return True, "Ring/chain system"

    return False, "Does not meet criteria for organic fundamental parent"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33245',
                          'name': 'organic fundamental parent',
                          'definition': 'An organic fundamental parent is a '
                                        'structure used as a basis for '
                                        'substitutive names in organic '
                                        'nomenclature, containing, in addition '
                                        'to one or more hydrogen atoms, a '
                                        'single atom of an element, a number '
                                        'of atoms (alike or different) linked '
                                        'together to form an unbranched chain, '
                                        'a monocyclic or polycyclic ring '
                                        'system, or a ring assembly or '
                                        'ring/chain system.',
                          'parents': ['CHEBI:37175', 'CHEBI:50860']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '[13:43:58] SMILES Parse Error: unclosed ring for input: '
             "'CC(C)C1=C/C=C(C)\\CC\\C=C(C)\\CC\\C=C(C)\\CC\x01'\n",
    'stdout': '',
    'num_true_positives': 188,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 1,
    'precision': 0.9038461538461539,
    'recall': 0.9947089947089947,
    'f1': 0.947103274559194,
    'accuracy': None}