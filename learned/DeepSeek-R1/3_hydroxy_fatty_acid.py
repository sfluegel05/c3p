"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI:134825 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxyl group (-OH) on the third carbon of the longest chain
    starting from the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find carboxylic acid groups (-COOH)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if not matches:
        return False, "No carboxylic acid group found"

    # Function to find the longest chain from the alpha carbon (carbon adjacent to COOH)
    def find_longest_chain(alpha_carbon):
        visited = set()
        queue = deque([(alpha_carbon, 0)])  # (atom, chain_length)
        max_length = 0
        while queue:
            atom, length = queue.popleft()
            if atom in visited:
                continue
            visited.add(atom)
            max_length = max(max_length, length)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != alpha_carbon.GetIdx():
                    queue.append((neighbor, length + 1))
        return max_length

    # Check each carboxylic acid group
    for match in matches:
        # The carbonyl carbon is match[0], the alpha carbon is connected to it
        carbonyl_carbon = mol.GetAtomWithIdx(match[0])
        alpha_carbons = [n for n in carbonyl_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if not alpha_carbons:
            continue  # No alpha carbon (unlikely)
        alpha_carbon = alpha_carbons[0]

        # Find the longest chain from alpha_carbon
        chain_length = find_longest_chain(alpha_carbon)
        if chain_length < 3:
            continue  # Chain too short to have a 3rd carbon

        # Traverse the chain to check for hydroxyl on the 3rd carbon
        current = alpha_carbon
        position = 1
        visited = set()
        while position < 3:
            visited.add(current.GetIdx())
            next_carbons = [n for n in current.GetNeighbors() 
                            if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
            if not next_carbons:
                break
            current = next_carbons[0]  # Follow the longest path (simplified)
            position += 1

        if position == 3:
            # Check if current (3rd carbon) has a hydroxyl group
            for neighbor in current.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    return True, "3-hydroxy group found on the main fatty acid chain"

    return False, "No 3-hydroxy group detected on the main chain from carboxylic acid"