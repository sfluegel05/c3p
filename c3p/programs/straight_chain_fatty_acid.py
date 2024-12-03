"""
Classifies: CHEBI:59202 straight-chain fatty acid
"""
from rdkit import Chem

def is_straight_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'O' in neighbors and neighbors.count('O') == 2:
                carboxyl_found = True
                break

    if not carboxyl_found:
        return False, "No carboxylic acid group found"

    # Check for unbranched carbon chain
    carbon_chain = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            if len([nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'C']) <= 2:
                carbon_chain.append(atom)

    if len(carbon_chain) < 2:
        return False, "No unbranched carbon chain found"

    # Check if the carbon chain is continuous and unbranched
    visited = set()
    def dfs(atom):
        visited.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == 'C' and nbr.GetIdx() not in visited:
                dfs(nbr)

    dfs(carbon_chain[0])
    if len(visited) != len(carbon_chain):
        return False, "Carbon chain is not continuous and unbranched"

    return True, "Straight-chain fatty acid"

# Example usage:
# print(is_straight_chain_fatty_acid("CCCCCCCCCCCCCCCCCC(O)=O"))  # Example SMILES string for testing


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59202',
                          'name': 'straight-chain fatty acid',
                          'definition': 'Any fatty acid whose skeletal carbon '
                                        'atoms form an unbranched open chain.',
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
    'stdout': '',
    'num_true_positives': 41,
    'num_false_positives': 16,
    'num_true_negatives': 4,
    'num_false_negatives': 2,
    'precision': 0.7192982456140351,
    'recall': 0.9534883720930233,
    'f1': 0.8200000000000001,
    'accuracy': None}