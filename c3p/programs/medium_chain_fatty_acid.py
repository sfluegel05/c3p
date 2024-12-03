"""
Classifies: CHEBI:59554 medium-chain fatty acid
"""
from rdkit import Chem

def is_medium_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acid (C6 to C12).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a carboxylic acid group
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check the chain length
    num_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            num_carbons += 1

    # Adjust the carbon count for carboxylic acid group
    if carboxylic_acid:
        num_carbons -= 1

    if 6 <= num_carbons <= 12:
        return True, f"Medium-chain fatty acid with {num_carbons} carbon atoms"
    else:
        return False, f"Chain length is {num_carbons} carbon atoms, not between 6 and 12"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59554',
                          'name': 'medium-chain fatty acid',
                          'definition': 'Any fatty acid with a chain length of '
                                        'between C6 and C12.',
                          'parents': ['CHEBI:35366']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '[01:30:26] SMILES Parse Error: syntax error while parsing: '
             'OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C=C\x01/[C@@H](C/C=C\\CC)O1)=O\n'
             '[01:30:26] SMILES Parse Error: Failed parsing SMILES '
             "'OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C=C\x01/[C@@H](C/C=C\\CC)O1)=O' "
             'for input: '
             "'OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C=C\x01/[C@@H](C/C=C\\CC)O1)=O'\n"
             '[01:30:26] SMILES Parse Error: syntax error while parsing: '
             'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O\n'
             '[01:30:26] SMILES Parse Error: Failed parsing SMILES '
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O' for input: "
             "'C(CC/C=C\\C=C\x01/O[C@H]1C/C=C\\CC)CCCCC(=O)O'\n",
    'stdout': '',
    'num_true_positives': 59,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 17,
    'precision': 1.0,
    'recall': 0.7763157894736842,
    'f1': 0.8740740740740741,
    'accuracy': None}