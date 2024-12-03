"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13 to C22).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

    # Check if the molecule is a fatty acid (contains a carboxyl group)
    carboxyl_group = Chem.MolFromSmarts('C(=O)O')
    has_carboxyl_group = mol.HasSubstructMatch(carboxyl_group)

    if not has_carboxyl_group:
        return False, "No carboxyl group found"

    # Check if the number of carbon atoms is within the range for long-chain fatty acids
    if 13 <= num_carbons <= 22:
        return True, f"Long-chain fatty acid with {num_carbons} carbon atoms"
    else:
        return False, f"Number of carbon atoms ({num_carbons}) is not in the range of C13 to C22"

# Example usage
smiles = "C(\CC)=C\C/C=C\C/C=C\C\C=C/C=C/C(CCCC(=O)O)O"  # 5-HEPE
print(is_long_chain_fatty_acid(smiles))  # Expected output: (True, "Long-chain fatty acid with 20 carbon atoms")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15904',
                          'name': 'long-chain fatty acid',
                          'definition': 'A fatty acid with a chain length '
                                        'ranging from C13 to C22.',
                          'parents': ['CHEBI:35366']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Long-chain fatty acid with 20 carbon atoms')\n",
    'num_true_positives': 179,
    'num_false_positives': 10,
    'num_true_negatives': 10,
    'num_false_negatives': 11,
    'precision': 0.9470899470899471,
    'recall': 0.9421052631578948,
    'f1': 0.9445910290237467,
    'accuracy': None}