"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide (a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a 2-furanone skeleton
    furanone_skeleton = Chem.MolFromSmarts('O=C1OC=CC1')
    if mol.HasSubstructMatch(furanone_skeleton):
        return True, "Contains 2-furanone skeleton"

    return False, "Does not contain 2-furanone skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50523',
                          'name': 'butenolide',
                          'definition': 'A gamma-lactone that consists of a '
                                        '2-furanone skeleton and its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:24129', 'CHEBI:37581']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:57:47] SMILES Parse Error: unclosed ring for input: '
             "'OC(=O)\\C=C1OC(=O)C=C\x01'\n"
             '[00:57:47] SMILES Parse Error: syntax error while parsing: '
             'O=C1O/C(=C\x02/[C@@H]([C@@H](CO)C)CC[C@H]2C)/C=C1C\n'
             '[00:57:47] SMILES Parse Error: Failed parsing SMILES '
             "'O=C1O/C(=C\x02/[C@@H]([C@@H](CO)C)CC[C@H]2C)/C=C1C' for input: "
             "'O=C1O/C(=C\x02/[C@@H]([C@@H](CO)C)CC[C@H]2C)/C=C1C'\n"
             '[00:57:47] SMILES Parse Error: syntax error while parsing: '
             'O1[C@@]([C@]2([C@@]\x03(C(O)(C(=O)[C@]([C@@H]2C(=O)/C=C/C=C/C)(C(=O)/C3=C(\\O)/C=C/C=C/C)C)C)[H])[H])(C(O)=C(C1=O)C)C\n'
             '[00:57:47] SMILES Parse Error: Failed parsing SMILES '
             "'O1[C@@]([C@]2([C@@]\x03(C(O)(C(=O)[C@]([C@@H]2C(=O)/C=C/C=C/C)(C(=O)/C3=C(\\O)/C=C/C=C/C)C)C)[H])[H])(C(O)=C(C1=O)C)C' "
             'for input: '
             "'O1[C@@]([C@]2([C@@]\x03(C(O)(C(=O)[C@]([C@@H]2C(=O)/C=C/C=C/C)(C(=O)/C3=C(\\O)/C=C/C=C/C)C)C)[H])[H])(C(O)=C(C1=O)C)C'\n"
             '[00:57:47] SMILES Parse Error: syntax error while parsing: '
             'CCC/C=C\x01/C=2CCCCC2C(=O)O1\n'
             '[00:57:47] SMILES Parse Error: Failed parsing SMILES '
             "'CCC/C=C\x01/C=2CCCCC2C(=O)O1' for input: "
             "'CCC/C=C\x01/C=2CCCCC2C(=O)O1'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 51,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}