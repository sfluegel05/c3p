"""
Classifies: CHEBI:50492 thiocarbonyl compound
"""
from rdkit import Chem

def is_thiocarbonyl_compound(smiles: str):
    """
    Determines if a molecule is a thiocarbonyl compound (contains the thiocarbonyl group, C=S).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiocarbonyl compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the thiocarbonyl group (C=S)
    thiocarbonyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'S' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                    thiocarbonyl_found = True
                    break
        if thiocarbonyl_found:
            break

    if thiocarbonyl_found:
        return True, "Thiocarbonyl group (C=S) found"
    else:
        return False, "Thiocarbonyl group (C=S) not found"

# Examples
print(is_thiocarbonyl_compound("COC1=CC=CC=C1NC(=S)NC2=CC=CC=N2"))  # Example 1
print(is_thiocarbonyl_compound("CC(C(C)(C)C)NC(=S)NC1=CC=CC(=C1)C(F)(F)F"))  # Example 2
print(is_thiocarbonyl_compound("CC1CCN(CC1)C(=S)NC2=CC=C(C=C2)F"))  # Example 3


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50492',
                          'name': 'thiocarbonyl compound',
                          'definition': 'Any compound containing the '
                                        'thiocarbonyl group, C=S.',
                          'parents': ['CHEBI:33261']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Thiocarbonyl group (C=S) found')\n"
              "(True, 'Thiocarbonyl group (C=S) found')\n"
              "(True, 'Thiocarbonyl group (C=S) found')\n",
    'num_true_positives': 30,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}