"""
Classifies: CHEBI:33249 organyl group
"""
from rdkit import Chem

def is_organyl_group(smiles: str):
    """
    Determines if a molecule is an organyl group (Any organic substituent group, regardless of functional type, having one free valence at a carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organyl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for one free valence at a carbon atom
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetNumImplicitHs() == 1:
            return True, "Molecule has one free valence at a carbon atom"
    
    return False, "Molecule does not have one free valence at a carbon atom"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33249',
                          'name': 'organyl group',
                          'definition': 'Any organic substituent group, '
                                        'regardless of functional type, having '
                                        'one free valence at a carbon atom.',
                          'parents': ['CHEBI:51447']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[22:58:07] SMILES Parse Error: syntax error while parsing: '
             'O(C[C@@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)CCCC(O)(C)C\n'
             '[22:58:07] SMILES Parse Error: Failed parsing SMILES '
             "'O(C[C@@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)CCCC(O)(C)C' "
             'for input: '
             "'O(C[C@@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)CCCC(O)(C)C'\n",
    'stdout': '',
    'num_true_positives': 14,
    'num_false_positives': 15,
    'num_true_negatives': 5,
    'num_false_negatives': 7,
    'precision': 0.4827586206896552,
    'recall': 0.6666666666666666,
    'f1': 0.56,
    'accuracy': None}