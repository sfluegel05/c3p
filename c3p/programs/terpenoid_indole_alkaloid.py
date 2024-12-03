"""
Classifies: CHEBI:65321 terpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_terpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a terpenoid indole alkaloid, defined as 'An indole alkaloid which is biosynthesised from L-tryptophan and isoprenoid building blocks.'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a terpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an indole ring
    indole_smarts = Chem.MolFromSmarts('C1=CC=C2C(=C1)C=CN2')
    if not mol.HasSubstructMatch(indole_smarts):
        return False, "No indole ring found"

    # Check for the presence of isoprenoid building blocks (C5 units)
    isoprenoid_smarts = Chem.MolFromSmarts('CC(C)=C')
    if not mol.HasSubstructMatch(isoprenoid_smarts):
        return False, "No isoprenoid building blocks found"

    # Check for the presence of L-tryptophan derived structure
    tryptophan_smarts = Chem.MolFromSmarts('C1=CC=C2C(=C1)C=CN2CC(C(=O)O)N')
    if not mol.HasSubstructMatch(tryptophan_smarts):
        return False, "No L-tryptophan derived structure found"

    return True, "Molecule is a terpenoid indole alkaloid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65321',
                          'name': 'terpenoid indole alkaloid',
                          'definition': 'An indole alkaloid which is '
                                        'biosynthesised from L-tryptophan and '
                                        'isoprenoid building blocks.',
                          'parents': ['CHEBI:38958']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[02:25:40] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@]4([C@@](C=C([C@]7(C[C@@]89C%10=CC=CC=C%10N%11C=C%12[C@]/%13(C[C@]%14([C@@]%15(CCN%14C\\C%13=C\\C)C%16=CC=CC=C%16N(C=C([C@]/%17(C[C@@]8(N7C\\C%17=C\\C)[H])[H])[C@@]9%11[H])[C@@]%12%15[H])[H])[H])[H])C6=O)([C@]1(C[C@@]32[H])[H])[H])[H]\n'
             '[02:25:40] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@]4([C@@](C=C([C@]7(C[C@@]89C%10=CC=CC=C%10N%11C=C%12[C@]/%13(C[C@]%14([C@@]%15(CCN%14C\\C%13=C\\C)C%16=CC=CC=C%16N(C=C([C@]/%17(C[C@@]8(N7C\\C%17=C\\C)[H])[H])[C@@]9%11[H])[C@@]%12%15[H])[H])[H])[H])C6=O)([C@]1(C[C@@]32[H])[H])[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@]4([C@@](C=C([C@]7(C[C@@]89C%10=CC=CC=C%10N%11C=C%12[C@]/%13(C[C@]%14([C@@]%15(CCN%14C\\C%13=C\\C)C%16=CC=CC=C%16N(C=C([C@]/%17(C[C@@]8(N7C\\C%17=C\\C)[H])[H])[C@@]9%11[H])[C@@]%12%15[H])[H])[H])[H])C6=O)([C@]1(C[C@@]32[H])[H])[H])[H]'\n"
             '[02:25:40] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2CC[C@@]3(C4=CC=CC=C4N5C=C([C@]1(C[C@]2([C@]35O)[H])[H])C(=O)OC)O\n'
             '[02:25:40] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2CC[C@@]3(C4=CC=CC=C4N5C=C([C@]1(C[C@]2([C@]35O)[H])[H])C(=O)OC)O' "
             'for input: '
             "'C/C=C\x01/CN2CC[C@@]3(C4=CC=CC=C4N5C=C([C@]1(C[C@]2([C@]35O)[H])[H])C(=O)OC)O'\n"
             '[02:25:40] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC=CC=C6N(C)[C@]35[H])[C@@]4(OC(\\C=C\\C7=CC(=C(C(=C7)OC)OC)OC)=O)[H])[H])C(=O)OC)[H])[H]\n'
             '[02:25:40] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC=CC=C6N(C)[C@]35[H])[C@@]4(OC(\\C=C\\C7=CC(=C(C(=C7)OC)OC)OC)=O)[H])[H])C(=O)OC)[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC=CC=C6N(C)[C@]35[H])[C@@]4(OC(\\C=C\\C7=CC(=C(C(=C7)OC)OC)OC)=O)[H])[H])C(=O)OC)[H])[H]'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 29,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}