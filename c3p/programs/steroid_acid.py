"""
Classifies: CHEBI:47891 steroid acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_steroid_acid(smiles: str):
    """
    Determines if a molecule is a steroid acid (any steroid substituted by at least one carboxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxyl group (COOH)
    carboxyl_group = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxyl_group):
        return False, "No carboxyl group found"

    # Check for the steroid structure (four fused rings: three 6-membered and one 5-membered)
    steroid_core = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCCC4')
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid core structure not found"

    return True, "Molecule is a steroid acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47891',
                          'name': 'steroid acid',
                          'definition': 'Any steroid substituted by at least '
                                        'one carboxy group.',
                          'parents': ['CHEBI:35341', 'CHEBI:64709']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:42:47] SMILES Parse Error: syntax error while parsing: '
             'O=C(CC[C@H](C=1[C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)C)C(O)(C)C\n'
             '[00:42:47] SMILES Parse Error: Failed parsing SMILES '
             "'O=C(CC[C@H](C=1[C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)C)C(O)(C)C' "
             'for input: '
             "'O=C(CC[C@H](C=1[C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)C)C(O)(C)C'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 51,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}