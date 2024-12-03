"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a phosphatidyl group
    phosphatidyl_smarts = '[O-]P(=O)(OCC[C@H](N)C(=O)O)OC[C@H](O[*])COC(=O)[*]'
    phosphatidyl_pattern = Chem.MolFromSmarts(phosphatidyl_smarts)
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "No phosphatidyl group esterified to serine detected"

    # Check for esterified hydroxy group of serine
    esterified_smarts = 'OC[C@H](O[*])COC(=O)[*]'
    esterified_pattern = Chem.MolFromSmarts(esterified_smarts)
    if not mol.HasSubstructMatch(esterified_pattern):
        return False, "No esterified hydroxy group of serine detected"

    return True, "Molecule is a phosphatidyl-L-serine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18303',
                          'name': 'phosphatidyl-L-serine',
                          'definition': 'A class of aminophospholipids in '
                                        'which a phosphatidyl group is '
                                        'esterified to the hydroxy group of '
                                        'serine.',
                          'parents': [   'CHEBI:52565',
                                         'CHEBI:60971',
                                         'CHEBI:84135']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:54:44] SMILES Parse Error: syntax error while parsing: '
             '[H][C@]12[C@H](C)\\C(O[C@]11CCC[N+]3([O-])CCC[C@@H](O1)[C@@]23[H])=C1OC(=O)C(C)=C\x01OC\n'
             '[19:54:44] SMILES Parse Error: Failed parsing SMILES '
             "'[H][C@]12[C@H](C)\\C(O[C@]11CCC[N+]3([O-])CCC[C@@H](O1)[C@@]23[H])=C1OC(=O)C(C)=C\x01OC' "
             'for input: '
             "'[H][C@]12[C@H](C)\\C(O[C@]11CCC[N+]3([O-])CCC[C@@H](O1)[C@@]23[H])=C1OC(=O)C(C)=C\x01OC'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 44,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}