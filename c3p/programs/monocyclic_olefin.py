"""
Classifies: CHEBI:36403 monocyclic olefin
"""
from rdkit import Chem

def is_monocyclic_olefin(smiles: str):
    """
    Determines if a molecule is a monocyclic olefin (a monocyclic hydrocarbon having any number of double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocyclic olefin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for exactly one ring
    if rings.NumRings() != 1:
        return False, "Not a monocyclic compound"

    # Check if the ring is a hydrocarbon (contains only C and H)
    ring_atoms = rings.AtomRings()[0]
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() not in {'C', 'H'}:
            return False, "Ring contains non-hydrocarbon atoms"

    # Check for the presence of double bonds in the ring
    has_double_bond = False
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_double_bond = True
                break

    if not has_double_bond:
        return False, "No double bonds in the ring"

    return True, "Monocyclic olefin"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36403',
                          'name': 'monocyclic olefin',
                          'definition': 'A monocyclic hydrocarbon having any '
                                        'number of double bonds.',
                          'parents': ['CHEBI:33642']},
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
    'num_true_positives': 10,
    'num_false_positives': 2,
    'num_true_negatives': 8,
    'num_false_negatives': 0,
    'precision': 0.8333333333333334,
    'recall': 1.0,
    'f1': 0.9090909090909091,
    'accuracy': None}