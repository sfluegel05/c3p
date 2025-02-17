"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: long-chain fatty acid (chain length from C13 to C22)
A long-chain fatty acid must have a carboxylic acid group and a linear carbon chain that
includes between 13 and 22 carbon atoms (counting the acid carbon as part of the chain).
Note: This is a heuristic method. In complex molecules or those having branching, the method
may overestimate the chain length.
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid must have a carboxylic acid group and possess an unbroken carbon chain 
    of length between 13 and 22 (inclusive), starting at the acid carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a carboxylic acid group.
    # We use a SMARTS pattern that matches the acid group: [CX3](=O)[OX2H]
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_group = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_group)
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # We'll use the acid carbon (first atom in the match) as the starting point of the chain.
    acid_idx = acid_matches[0][0]

    # Build a set of atom indices that correspond to carbon atoms.
    carbon_indices = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    if acid_idx not in carbon_indices:
        # In principle the acid carbon should be carbon.
        return False, "Acid group does not contain a carbon atom as expected"

    # Build an adjacency list for the carbon subgraph.
    # We iterate over each carbon atom and record carbon neighbors
    carbon_neighbors = {idx: [] for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors[idx].append(neighbor.GetIdx())

    # Now find the longest simple path (chain) of carbons starting from the acid carbon.
    # DFS-based recursion to compute maximum chain length (number of carbon atoms in path)
    def dfs(current, visited):
        max_length = 1  # count current atom
        for nbr in carbon_neighbors.get(current, []):
            if nbr not in visited:
                candidate = 1 + dfs(nbr, visited | {nbr})
                if candidate > max_length:
                    max_length = candidate
        return max_length

    longest_chain_length = dfs(acid_idx, {acid_idx})

    # Check if the longest chain length is between 13 and 22 carbons (inclusive)
    if longest_chain_length < 13:
        return False, f"Longest carbon chain is {longest_chain_length} carbons, which is too short"
    if longest_chain_length > 22:
        return False, f"Longest carbon chain is {longest_chain_length} carbons, which is too long"

    return True, f"Contains a carboxylic acid group with a connected carbon chain of {longest_chain_length} carbons"

# Example (uncomment to test):
# test_smiles = "OC(=O)CCCCCCCCCCCCCCC([2H])([2H])[2H]"  # palmitic acid(d3)
# result, reason = is_long_chain_fatty_acid(test_smiles)
# print(result, reason)