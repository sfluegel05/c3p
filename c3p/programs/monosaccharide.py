"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carbonyl group (aldehyde or ketone)
    carbonyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and any([neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0 for neighbor in atom.GetNeighbors()]):
            carbonyl_found = True
            break

    if not carbonyl_found:
        return False, "No carbonyl group found"

    # Check for the presence of at least three carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count < 3:
        return False, "Less than three carbon atoms found"

    # Check for polyhydroxy groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any([neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()]))
    if hydroxyl_count < 2:
        return False, "Less than two hydroxyl groups found"

    # Check for glycosidic connections
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) > 1:
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O' and any([neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in ring for neighbor in atom.GetNeighbors()]):
                    return False, "Glycosidic connections to other units found"

    return True, "Monosaccharide structure identified"

# Example usage
smiles = "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"
result, reason = is_monosaccharide(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35381',
                          'name': 'monosaccharide',
                          'definition': 'Parent monosaccharides are '
                                        'polyhydroxy aldehydes '
                                        'H[CH(OH)]nC(=O)H or polyhydroxy '
                                        'ketones H-[CHOH]n-C(=O)[CHOH]m-H with '
                                        'three or more carbon atoms. The '
                                        "generic term 'monosaccharide' (as "
                                        'opposed to oligosaccharide or '
                                        'polysaccharide) denotes a single '
                                        'unit, without glycosidic connection '
                                        'to other such units. It includes '
                                        'aldoses, dialdoses, aldoketoses, '
                                        'ketoses and diketoses, as well as '
                                        'deoxy sugars, provided that the '
                                        'parent compound has a (potential) '
                                        'carbonyl group.',
                          'parents': ['CHEBI:16646']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'False No carbonyl group found\n',
    'num_true_positives': 33,
    'num_false_positives': 13,
    'num_true_negatives': 7,
    'num_false_negatives': 24,
    'precision': 0.717391304347826,
    'recall': 0.5789473684210527,
    'f1': 0.6407766990291262,
    'accuracy': None}