"""
Classifies: CHEBI:38092 cardenolide glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cardenolide_glycoside(smiles: str):
    """
    Determines if a molecule is a cardenolide glycoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardenolide glycoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the cardenolide core structure (steroid backbone with lactone ring)
    cardenolide_core_smarts = '[C@]12CC[C@]3([H])[C@@]1(CC[C@]4([H])[C@@]3(CC[C@@H](C4)O)C)C'
    cardenolide_core = Chem.MolFromSmarts(cardenolide_core_smarts)
    
    if not mol.HasSubstructMatch(cardenolide_core):
        return False, "Molecule does not contain cardenolide core structure"

    # Check for glycosylation at position 3
    glycoside_smarts = '[C@]1([H])[C@H](O)[C@@H](O[C@H]2[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O2)C)[C@H](C)C[C@@H]1O'
    glycoside = Chem.MolFromSmarts(glycoside_smarts)
    
    if not mol.HasSubstructMatch(glycoside):
        return False, "Molecule does not have glycosylation at position 3"

    return True, "Molecule is a cardenolide glycoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38092',
                          'name': 'cardenolide glycoside',
                          'definition': 'Any member of the class of '
                                        'cardenolides with glycosyl residues '
                                        'attached to position 3.',
                          'parents': ['CHEBI:74634', 'CHEBI:83970']},
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