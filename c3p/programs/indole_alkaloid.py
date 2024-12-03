"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid (an alkaloid containing an indole skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the indole substructure
    indole_smarts = 'c1ccc2c(c1)[nH]c2'
    indole_mol = Chem.MolFromSmarts(indole_smarts)

    # Check if the molecule contains the indole substructure
    if mol.HasSubstructMatch(indole_mol):
        return True, "Molecule contains an indole skeleton"
    else:
        return False, "Molecule does not contain an indole skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38958',
                          'name': 'indole alkaloid',
                          'definition': 'An alkaloid containing an indole '
                                        'skeleton.',
                          'parents': ['CHEBI:22315']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]' "
             'for input: '
             "'C/C=C\x01/C2CC[N@+]3(CC[C@]4(C5=CC=CC=C5NC4=C2C(=O)OC)[C@@]13[H])[O-]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2[C@]3(CC=4C5=CC=CC=C5NC4[C@@]2(C[C@@]1([C@H]3C(=O)OC)[H])[H])[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2[C@]3(CC=4C5=CC=CC=C5NC4[C@@]2(C[C@@]1([C@H]3C(=O)OC)[H])[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN2[C@]3(CC=4C5=CC=CC=C5NC4[C@@]2(C[C@@]1([C@H]3C(=O)OC)[H])[H])[H]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@]36OC(C7=C4C=CC(=C7O)O)=O)[H])[H])(C(=O)OC)[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@]36OC(C7=C4C=CC(=C7O)O)=O)[H])[H])(C(=O)OC)[H]' "
             'for input: '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@]36OC(C7=C4C=CC(=C7O)O)=O)[H])[H])(C(=O)OC)[H]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC(=CC=C6N(C)[C@]35[H])O)[C@@]4(OC(\\C=C\\C7=CC(=C(C(=C7)OC)OC)OC)=O)[H])[H])C(=O)OC)[H])[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC(=CC=C6N(C)[C@]35[H])O)[C@@]4(OC(\\C=C\\C7=CC(=C(C(=C7)OC)OC)OC)=O)[H])[H])C(=O)OC)[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC(=CC=C6N(C)[C@]35[H])O)[C@@]4(OC(\\C=C\\C7=CC(=C(C(=C7)OC)OC)OC)=O)[H])[H])C(=O)OC)[H])[H]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC=CC=C6N(C)[C@]35[H])[C@@H]4O)[H])C(=O)OC)[H])[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC=CC=C6N(C)[C@]35[H])[C@@H]4O)[H])C(=O)OC)[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC=CC=C6N(C)[C@]35[H])[C@@H]4O)[H])C(=O)OC)[H])[H]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@@]36OC7=C4C=C8C(=C7O)C(=O)O[C@@]9%10[C@@]%11(C[C@]%12(/C(=C\\C)/CN%11CC[C@@]8%10C%13=CC=CC=C%13N9[C@@]%12(C(=O)OC)[H])[H])[H])[H])[H])(C(=O)OC)[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@@]36OC7=C4C=C8C(=C7O)C(=O)O[C@@]9%10[C@@]%11(C[C@]%12(/C(=C\\C)/CN%11CC[C@@]8%10C%13=CC=CC=C%13N9[C@@]%12(C(=O)OC)[H])[H])[H])[H])[H])(C(=O)OC)[H]' "
             'for input: '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@@]36OC7=C4C=C8C(=C7O)C(=O)O[C@@]9%10[C@@]%11(C[C@]%12(/C(=C\\C)/CN%11CC[C@@]8%10C%13=CC=CC=C%13N9[C@@]%12(C(=O)OC)[H])[H])[H])[H])[H])(C(=O)OC)[H]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/C[N+]2(C)CCC=3C4=CC=CC=C4N5C(CO)[C@]1(C[C@]2(C35)[H])[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/C[N+]2(C)CCC=3C4=CC=CC=C4N5C(CO)[C@]1(C[C@]2(C35)[H])[H]' "
             'for input: '
             "'C/C=C\x01/C[N+]2(C)CCC=3C4=CC=CC=C4N5C(CO)[C@]1(C[C@]2(C35)[H])[H]'\n"
             '[00:25:30] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN(C)CC[C@@]23C4=CC=CC=C4N5[C@]3([C@](CO[C@@H](C5=O)O)([C@]1(CC2=O)[H])[H])[H]\n'
             '[00:25:30] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN(C)CC[C@@]23C4=CC=CC=C4N5[C@]3([C@](CO[C@@H](C5=O)O)([C@]1(CC2=O)[H])[H])[H]' "
             'for input: '
             "'C/C=C\x01/CN(C)CC[C@@]23C4=CC=CC=C4N5[C@]3([C@](CO[C@@H](C5=O)O)([C@]1(CC2=O)[H])[H])[H]'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 51,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}