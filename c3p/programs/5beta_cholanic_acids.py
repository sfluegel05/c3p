"""
Classifies: CHEBI:36248 5beta-cholanic acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_5beta_cholanic_acids(smiles: str):
    """
    Determines if a molecule is a 5beta-cholanic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 5beta-cholanic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the 5beta-cholane skeleton
    cholane_skeleton = Chem.MolFromSmarts('C1[C@H](C[C@@]2([C@](C1)([C@@]3([C@@]([C@H](CC2)C)(CC4)[H])C)[H])[H])[H]')
    if not mol.HasSubstructMatch(cholane_skeleton):
        return False, "No 5beta-cholane skeleton found"

    # Check for carboxylic acid functionality
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid functionality found"

    return True, "Molecule is a 5beta-cholanic acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36248',
                          'name': '5beta-cholanic acids',
                          'definition': 'Members of the class of cholanic '
                                        'acids based on a 5beta-cholane '
                                        'skeleton.',
                          'parents': ['CHEBI:136889', 'CHEBI:36278']},
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