"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 3beta-hydroxy group
    pattern_3beta_hydroxy = Chem.MolFromSmarts("[C@H](O)[C@@H]1CC[C@H]2")
    if not mol.HasSubstructMatch(pattern_3beta_hydroxy):
        return False, "No 3beta-hydroxy group found"

    # Check for double bond between positions 5 and 6
    pattern_Delta_5 = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(pattern_Delta_5)
    double_bond_positions = [(match[0], match[1]) for match in double_bond_matches]
    
    # Assuming positions 5 and 6 are part of a steroid framework
    steroid_positions = [(4, 5), (5, 6), (6, 7), (7, 8)]
    if not any(pos in steroid_positions for pos in double_bond_positions):
        return False, "No double bond between positions 5 and 6 found"
    
    return True, "Molecule is a 3beta-hydroxy-Delta(5)-steroid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:1722',
                          'name': '3beta-hydroxy-Delta(5)-steroid',
                          'definition': 'Any 3beta-hydroxy-steroid that '
                                        'contains a double bond between '
                                        'positions 5 and 6.',
                          'parents': ['CHEBI:36836']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': 'Python argument types in\n'
             '    Mol.HasSubstructMatch(Mol, NoneType)\n'
             'did not match C++ signature:\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'RDKit::SubstructMatchParameters params=True)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'RDKit::SubstructMatchParameters params)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}