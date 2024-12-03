"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define functional groups for conjugation
    conjugation_groups = [
        'C(C(=O)O)N',  # Glycine
        'C(C(=O)O)NC(CS(=O)(=O)O)',  # Taurine
        'S(=O)(=O)O',  # Sulfuric acid (sulfate)
        'C(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O',  # Glucuronic acid (glucuronate)
        'C(CO)O',  # Glucose
        'CC(=O)SCCNC(=O)C',  # Coenzyme A
    ]

    # Check for bile acid core structure
    bile_acid_core = Chem.MolFromSmarts('C1[C@@]2([C@]3([C@@]([C@]4([C@@]([C@@]5([C@](C[C@@H]3O)(C[C@H](O)CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H]')
    if not mol.HasSubstructMatch(bile_acid_core):
        return False, "No bile acid core structure found"

    # Check for conjugation groups
    for group in conjugation_groups:
        group_mol = Chem.MolFromSmarts(group)
        if mol.HasSubstructMatch(group_mol):
            return True, f"Bile acid conjugate with group: {group}"

    return False, "No conjugation groups found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36249',
                          'name': 'bile acid conjugate',
                          'definition': 'Any bile acid conjugated to a '
                                        'functional group that gives '
                                        'additional hydrophilicity or charge '
                                        'to the molecule. Molecules used for '
                                        'conjugation are: glycine, taurine '
                                        '(and other amino acids); sulfuric '
                                        "acid (for which the term ''sulfate'' "
                                        'may be used); glucuronic acid (for '
                                        "which the term ''glucuronate'' may be "
                                        'used); glucose and other uncharged '
                                        'sugars; and coenzyme A.',
                          'parents': ['CHEBI:36078']},
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