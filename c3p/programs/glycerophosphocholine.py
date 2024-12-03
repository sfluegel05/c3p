"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for choline group
    choline_pattern = Chem.MolFromSmarts("N+(C)(C)C")
    if not mol.HasSubstructMatch(choline_pattern):
        return False, "No choline group found"

    return True, "Molecule is a glycerophosphocholine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36313',
                          'name': 'glycerophosphocholine',
                          'definition': 'The glycerol phosphate ester of a '
                                        'phosphocholine. A nutrient with many '
                                        'different roles in human health.',
                          'parents': ['CHEBI:36700', 'CHEBI:37739']},
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