"""
Classifies: CHEBI:35568 mancude ring
"""
from rdkit import Chem

def is_mancude_ring(smiles: str):
    """
    Determines if a molecule is a mancude ring (a ring having the maximum number of noncumulative double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mancude ring, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check if there is at least one ring
    if not rings.AtomRings():
        return False, "No rings found"

    # Check for mancude rings
    for ring in rings.BondRings():
        # Check for maximum number of noncumulative double bonds
        double_bonds = 0
        single_bonds = 0
        for bond_idx in ring:
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bonds += 1
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                single_bonds += 1

        # Mancude ring condition: the number of double bonds should be equal to the number of single bonds
        if double_bonds == single_bonds:
            return True, "Mancude ring found"

    return False, "No mancude ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35568',
                          'name': 'mancude ring',
                          'definition': 'Any molecular entity that consists of '
                                        'a ring having (formally) the maximum '
                                        'number of noncumulative double bonds.',
                          'parents': ['CHEBI:23367']},
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
    'num_true_positives': 12,
    'num_false_positives': 10,
    'num_true_negatives': 5,
    'num_false_negatives': 3,
    'precision': 0.5454545454545454,
    'recall': 0.8,
    'f1': 0.6486486486486486,
    'accuracy': None}