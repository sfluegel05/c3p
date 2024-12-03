"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline (dihydropyrrole).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all 5-membered rings containing nitrogen
    pyrroline_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'N' for atom in atoms):
                pyrroline_rings.append(ring)

    if not pyrroline_rings:
        return False, "No 5-membered rings containing nitrogen found"

    # Check if the ring is partially saturated (dihydropyrrole)
    for ring in pyrroline_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring and bond.GetBondType() == Chem.BondType.DOUBLE)
        single_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring and bond.GetBondType() == Chem.BondType.SINGLE)
        if double_bonds == 1 and single_bonds == 4:
            return True, "Dihydropyrrole ring found"
        elif double_bonds == 2 and single_bonds == 3:
            return True, "Dihydropyrrole ring found"
        # Additional check for partially saturated rings
        elif double_bonds == 0 and single_bonds == 5:
            return True, "Dihydropyrrole ring found"

    return False, "No dihydropyrrole ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23763',
                          'name': 'pyrroline',
                          'definition': 'Any organic heteromonocyclic compound '
                                        'with a structure based on a '
                                        'dihydropyrrole.',
                          'parents': ['CHEBI:25693', 'CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '[20:21:08] SMILES Parse Error: syntax error while parsing: '
             'CC1=CC=C(/C=C\x02/CCN=C2)O1\n'
             '[20:21:08] SMILES Parse Error: Failed parsing SMILES '
             "'CC1=CC=C(/C=C\x02/CCN=C2)O1' for input: "
             "'CC1=CC=C(/C=C\x02/CCN=C2)O1'\n",
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 1,
    'num_true_negatives': 11,
    'num_false_negatives': 1,
    'precision': 0.9166666666666666,
    'recall': 0.9166666666666666,
    'f1': 0.9166666666666666,
    'accuracy': None}