"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an indole structure
    indole_smarts = 'c1c[nH]c2ccccc12'
    indole = Chem.MolFromSmarts(indole_smarts)
    if not mol.HasSubstructMatch(indole):
        return False, "No indole structure found"

    # Check for terpenoid structure (isoprenoid units)
    isoprenoid_smarts = 'C=C(C)C'
    isoprenoid = Chem.MolFromSmarts(isoprenoid_smarts)
    if not mol.HasSubstructMatch(isoprenoid):
        return False, "No isoprenoid units found"

    return True, "Molecule is a monoterpenoid indole alkaloid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65323',
                          'name': 'monoterpenoid indole alkaloid',
                          'definition': 'A terpenoid indole alkaloid which is '
                                        'biosynthesised from L-tryptophan and '
                                        'diisoprenoid (usually secolaganin) '
                                        'building blocks.',
                          'parents': ['CHEBI:65321']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '[02:26:26] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]\n'
             '[02:26:26] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]' "
             'for input: '
             "'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]'\n"
             '[02:26:26] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]\n'
             '[02:26:26] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]'\n"
             '[02:26:26] SMILES Parse Error: syntax error while parsing: '
             'C(=C\x01/CN2CCC3=C([C@@]2(C[C@@]1(CC=4C5=C(C=6C=CC=CC6N5)C=CN4)[H])[H])NC=7C=CC(=CC37)O)\\C\n'
             '[02:26:26] SMILES Parse Error: Failed parsing SMILES '
             "'C(=C\x01/CN2CCC3=C([C@@]2(C[C@@]1(CC=4C5=C(C=6C=CC=CC6N5)C=CN4)[H])[H])NC=7C=CC(=CC37)O)\\C' "
             'for input: '
             "'C(=C\x01/CN2CCC3=C([C@@]2(C[C@@]1(CC=4C5=C(C=6C=CC=CC6N5)C=CN4)[H])[H])NC=7C=CC(=CC37)O)\\C'\n"
             '[02:26:26] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2CCC=3C=4C=C(C=CC4NC3C2CC1CCO)O\n'
             '[02:26:26] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2CCC=3C=4C=C(C=CC4NC3C2CC1CCO)O' for input: "
             "'C/C=C\x01/CN2CCC=3C=4C=C(C=CC4NC3C2CC1CCO)O'\n",
    'stdout': '',
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 22,
    'precision': 1.0,
    'recall': 0.043478260869565216,
    'f1': 0.08333333333333333,
    'accuracy': None}