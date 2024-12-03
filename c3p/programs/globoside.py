"""
Classifies: CHEBI:61360 globoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_globoside(smiles: str):
    """
    Determines if a molecule is a globoside (glycosphingolipid with an N-acetylgalactosaminyl residue at the non-reducing end).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a globoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for N-acetylgalactosaminyl residue
    pattern = Chem.MolFromSmarts("NC(C)=O")
    if not mol.HasSubstructMatch(pattern):
        return False, "No N-acetylgalactosaminyl residue found"

    # Check for glycosphingolipid structure
    glycosphingolipid_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O[C@H]2[C@@H](O)[C@H](O[C@H]3[C@@H](O)[C@H](O[C@H]4[C@@H](O)[C@H](O[C@H]5[C@@H](O)[C@H](O[C@H]6[C@@H](O)[C@H](O)[C@H](O6)CO)O)O)O)O)O[C@@H]5O[C@H](CO)[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycosphingolipid_pattern):
        return False, "No glycosphingolipid structure found"

    return True, "Molecule is a globoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61360',
                          'name': 'globoside',
                          'definition': 'Any glycosphingolipid where the '
                                        'oligosaccharide component has an '
                                        'N-acetylgalactosaminyl residue at the '
                                        'non-reducing end.',
                          'parents': ['CHEBI:24402']},
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