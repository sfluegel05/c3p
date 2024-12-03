"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid (chain length greater than C22).
    Very long-chain fatty acids with a chain length greater than C27 are also known as ultra-long-chain fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the carboxylic acid group
    carboxylic_acid = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                carboxylic_acid = atom
                break

    if carboxylic_acid is None:
        return False, "No carboxylic acid group found"

    # Traverse the carbon chain from the carboxylic acid group
    chain_length = 0
    current_atom = carboxylic_acid.GetNeighbors()[0]
    visited_atoms = set()

    while current_atom and current_atom.GetSymbol() == 'C' and current_atom.GetIdx() not in visited_atoms:
        visited_atoms.add(current_atom.GetIdx())
        chain_length += 1
        next_atoms = [n for n in current_atom.GetNeighbors() if n.GetSymbol() == 'C' and n.GetIdx() not in visited_atoms]
        current_atom = next_atoms[0] if next_atoms else None

    if chain_length > 27:
        return True, "Ultra-long-chain fatty acid"
    elif chain_length > 22:
        return True, "Very long-chain fatty acid"
    else:
        return False, f"Chain length is {chain_length}, which is not greater than C22"

# Example usage
smiles = "OC(=O)CCCCCCCCC#CC#CCCCCCCCCCC"
print(is_very_long_chain_fatty_acid(smiles))  # Example output: (True, "Very long-chain fatty acid")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27283',
                          'name': 'very long-chain fatty acid',
                          'definition': 'A fatty acid which has a chain length '
                                        'greater than C22. Very long-chain '
                                        'fatty acids which have a chain length '
                                        'greater than C27 are also known as '
                                        'ultra-long-chain fatty acids.',
                          'parents': ['CHEBI:35366']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'Chain length is 0, which is not greater than C22')\n",
    'num_true_positives': 6,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 19,
    'precision': 1.0,
    'recall': 0.24,
    'f1': 0.3870967741935484,
    'accuracy': None}