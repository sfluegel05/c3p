"""
Classifies: CHEBI:60926 amino monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_amino_monosaccharide(smiles: str):
    """
    Determines if a molecule is an amino monosaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino monosaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a monosaccharide
    # The molecule should contain a ring with at least 5 or 6 members and oxygen atoms
    ring_info = mol.GetRingInfo()
    ring_found = False
    for ring in ring_info.AtomRings():
        if 5 <= len(ring) <= 6:
            if any(mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'O' for atom_idx in ring):
                ring_found = True
                break

    if not ring_found:
        return False, "No monosaccharide ring found"

    # Check for the presence of an amino group replacing a hydroxy group
    amino_group_present = False
    hydroxy_group_replaced = False

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for n_neighbor in neighbor.GetNeighbors():
                        if n_neighbor.GetSymbol() == 'O':
                            amino_group_present = True
                            break

    if amino_group_present:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'O':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C':
                        for n_neighbor in neighbor.GetNeighbors():
                            if n_neighbor.GetSymbol() == 'N':
                                hydroxy_group_replaced = True
                                break

    if amino_group_present and hydroxy_group_replaced:
        return True, "Amino monosaccharide with one hydroxy group replaced by an amino group"

    return False, "Does not meet amino monosaccharide criteria"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60926',
                          'name': 'amino monosaccharide',
                          'definition': 'Any  amino sugar that is a '
                                        'monosaccharide in which one alcoholic '
                                        'hydroxy group is replaced by an amino '
                                        'group.',
                          'parents': ['CHEBI:28963']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[01:44:18] SMILES Parse Error: syntax error while parsing: '
             'O=C(N[C@H]1[C@H](O[C@H](CO)[C@H]([C@@H]1O)O)OC/C(=C\x02/[C@@H]3[C@H]([C@@H](C)CC3)CC(C2)(C)C)/CO)C\n'
             '[01:44:18] SMILES Parse Error: Failed parsing SMILES '
             "'O=C(N[C@H]1[C@H](O[C@H](CO)[C@H]([C@@H]1O)O)OC/C(=C\x02/[C@@H]3[C@H]([C@@H](C)CC3)CC(C2)(C)C)/CO)C' "
             'for input: '
             "'O=C(N[C@H]1[C@H](O[C@H](CO)[C@H]([C@@H]1O)O)OC/C(=C\x02/[C@@H]3[C@H]([C@@H](C)CC3)CC(C2)(C)C)/CO)C'\n",
    'stdout': '',
    'num_true_positives': 22,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 5,
    'precision': 0.5238095238095238,
    'recall': 0.8148148148148148,
    'f1': 0.6376811594202898,
    'accuracy': None}