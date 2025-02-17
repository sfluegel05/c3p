"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated Fatty Acid

Definition: A fatty acid must have a terminal carboxylic acid group and an acyclic carbon chain 
(with at least 5 carbons overall). A monounsaturated fatty acid (MUFA) has exactly one carbon–carbon 
double bond (and no triple bonds) along its longest continuous acyl chain.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES string.
      2. Locates a carboxylic acid group (accepting both protonated and deprotonated forms) by SMARTS.
      3. Requires that the acid carbon is terminal, i.e. it has exactly one carbon neighbour.
      4. From that carbon neighbour (chain_start), constructs a graph of carbon atoms connected by 
         non‐ring bonds (thus following the acyclic fatty acyl chain).
      5. Uses DFS to extract the longest continuous chain starting from chain_start.
      6. Checks that the extracted chain is long enough (at least 4 carbons, so that with the acid carbon
         there are at least 5 carbons total).
      7. Walks along the consecutive bonds in the chain counting double and triple bonds.
      8. Returns True if exactly one double bond (and no triple bonds) is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of True and a success explanation if the molecule qualifies as a
                     monounsaturated fatty acid; otherwise False with an explanation.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Locate a carboxylic acid group.
    # The SMARTS here matches both protonated and deprotonated acids.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Molecule does not contain a carboxylic acid group, thus not a fatty acid"
    
    # Assume the first match is our target acid group.
    acid_match = acid_matches[0]
    acid_idx = acid_match[0]  # acid carbon: first atom in pattern.
    acid_atom = mol.GetAtomWithIdx(acid_idx)
    
    # The acid carbon must be terminal: it should have exactly one carbon neighbour.
    carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Acid group is not terminal; acid carbon does not have exactly one carbon neighbour"
    
    chain_start_atom = carbon_neighbors[0]
    
    # Check overall that there are enough carbons in the molecule (at least 5 overall).
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 5:
        return False, "Not enough carbon atoms to be considered a fatty acid"
    
    # Build a graph of non‐ring bonds connecting carbon atoms.
    # We restrict to bonds that are NOT in any ring.
    graph = {}
    # Get all carbon atom indices.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    for idx in carbon_indices:
        graph[idx] = []
    # Add edges for non‐ring C–C bonds.
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if a.GetAtomicNum() == 6 and b.GetAtomicNum() == 6:
            if not bond.IsInRing():
                graph[a.GetIdx()].append(b.GetIdx())
                graph[b.GetIdx()].append(a.GetIdx())
    
    # Use DFS to find the longest path (chain) starting from the chain_start_atom.
    longest_path = []
    def dfs(current, path, visited):
        nonlocal longest_path
        # Update longest_path if we found a longer path.
        if len(path) > len(longest_path):
            longest_path = path.copy()
        for neighbor in graph.get(current, []):
            # We only follow the chain in a “linear” direction.
            if neighbor not in visited:
                visited.add(neighbor)
                path.append(neighbor)
                dfs(neighbor, path, visited)
                path.pop()
                visited.remove(neighbor)
    
    start_idx = chain_start_atom.GetIdx()
    # Ensure the starting atom is part of our graph.
    if start_idx not in graph:
        return False, "Chain-start atom is not present in the carbon graph"
    dfs(start_idx, [start_idx], {start_idx})
    
    # For a fatty acid chain (excluding the acid carbon), we require a minimum length.
    # (The overall molecule must have at least 5 carbons: acid carbon + chain of at least 4.)
    if len(longest_path) < 4:
        return False, "Fatty acid chain too short; expected at least 4 carbons in the chain"
    
    # Count unsaturations (double and triple bonds) along the consecutive bonds in the longest chain.
    double_count = 0
    triple_count = 0
    for i in range(len(longest_path) - 1):
        bond = mol.GetBondBetweenAtoms(longest_path[i], longest_path[i+1])
        if bond is None:
            continue
        bt = bond.GetBondType()
        if bt == Chem.BondType.DOUBLE:
            double_count += 1
        elif bt == Chem.BondType.TRIPLE:
            triple_count += 1
    
    if triple_count > 0:
        return False, "Fatty acid chain contains a triple bond, which is not allowed"
    
    if double_count != 1:
        return False, f"Fatty acid chain unsaturation count is {double_count} instead of exactly 1"
    
    return True, ("Molecule contains a terminal carboxylic acid group and an acyl carbon chain "
                  "with exactly one unsaturation (double bond) and no triple bonds, classifying it as a monounsaturated fatty acid")

# Example usage:
# Uncomment the following lines to test with one of the provided SMILES:
# if __name__ == '__main__':
#     test_smiles = "CCCCCCCC\\C=C/CCCCCCC(O)=O"  # cis-8-heptadecenoic acid
#     result, reason = is_monounsaturated_fatty_acid(test_smiles)
#     print(result, reason)