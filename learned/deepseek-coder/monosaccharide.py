"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or ketone with at least three carbon atoms,
    and it should not be part of a larger oligosaccharide or polysaccharide structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least three carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Less than three carbon atoms, not a monosaccharide"

    # Check for polyhydroxy structure (at least two hydroxyl groups)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 2:
        return False, "Insufficient hydroxyl groups for a monosaccharide"

    # Check for a carbonyl group (aldehyde or ketone)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group (aldehyde or ketone) found"

    # Check for glycosidic bonds (indicating oligo/polysaccharide)
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][OX2]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Glycosidic bond detected, likely part of an oligo/polysaccharide"

    # Check for ring structure (common in monosaccharides)
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No ring structure found, uncommon for monosaccharides"

    # Check for typical monosaccharide ring size (5 or 6 atoms)
    ring_sizes = set()
    for ring in ring_info.AtomRings():
        ring_sizes.add(len(ring))
    if not (5 in ring_sizes or 6 in ring_sizes):
        return False, "Ring size not typical for monosaccharides (expected 5 or 6 atoms)"

    return True, "Polyhydroxy aldehyde or ketone with at least three carbon atoms, no glycosidic bonds, and typical ring structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35381',
                          'name': 'monosaccharide',
                          'definition': 'Parent monosaccharides are polyhydroxy '
                                        'aldehydes H[CH(OH)]nC(=O)H or '
                                        'polyhydroxy ketones H-[CHOH]n-C(=O)[CHOH]m-H '
                                        'with three or more carbon atoms. The generic '
                                        'term \'monosaccharide\' (as opposed to '
                                        'oligosaccharide or polysaccharide) denotes a '
                                        'single unit, without glycosidic connection to '
                                        'other such units. It includes aldoses, '
                                        'dialdoses, aldoketoses, ketoses and diketoses, '
                                        'as well as deoxy sugars, provided that the '
                                        'parent compound has a (potential) carbonyl '
                                        'group.',
                          'parents': ['CHEBI:33838', 'CHEBI:33839']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}