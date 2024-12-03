"""
Classifies: CHEBI:52396 alpha-oxyketone
"""
from rdkit import Chem

def is_alpha_oxyketone(smiles: str):
    """
    Determines if a molecule is an alpha-oxyketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-oxyketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carbonyl group (C=O)
    carbonyl_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2:
                    carbonyl_atoms.append(atom)
                    break

    if not carbonyl_atoms:
        return False, "No carbonyl group found"

    # Check if the carbonyl carbon is bonded to an oxygen atom (oxy group) and ensure it's not bonded to hydrogen
    for carbonyl_atom in carbonyl_atoms:
        oxy_group = False
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 1:
                oxy_group = True
            if neighbor.GetSymbol() == 'H':
                return False, "Carbonyl carbon is bonded to hydrogen atoms"

        if oxy_group:
            return True, "Molecule is an alpha-oxyketone"

    return False, "No oxy group found bonded to the carbonyl carbon"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52396',
                          'name': 'alpha-oxyketone',
                          'definition': 'An oxyketone with the general formula '
                                        'R2C(=O) (R=/=H) where one or more of '
                                        'the R groups contains an oxy (-O-) '
                                        'group and the oxy and carbonyl groups '
                                        'are bonded to the same carbon atom.',
                          'parents': ['CHEBI:52395']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '[01:15:39] SMILES Parse Error: unclosed ring for input: '
             "'C[C@H]1C\\C=C\\[C@H]2[C@H](O)C(C)=C(C)[C@H]3[C@H](Cc4c[nH]c5ccccc45)NC(=O)[C@@]23C(=O)CC[C@H](O)C(=O)\\C(C)=C\x01'\n"
             '[01:15:39] SMILES Parse Error: unclosed ring for input: '
             "'C[C@H]1C\\C=C\\[C@H]2[C@H](O)C(C)=C(C)[C@H]3[C@H](Cc4c[nH]c5ccccc45)NC(=O)[C@@]23C(=O)\\C=C\\C(=O)[C@H](O)\\C(C)=C\x01'\n",
    'stdout': '',
    'num_true_positives': 17,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 34,
    'precision': 0.6071428571428571,
    'recall': 0.3333333333333333,
    'f1': 0.430379746835443,
    'accuracy': None}