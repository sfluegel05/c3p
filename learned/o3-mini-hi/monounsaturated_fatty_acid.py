"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated Fatty Acid

Definition: A fatty acid is defined to have a terminal carboxylic acid group and an acyclic carbon chain 
(with at least 5 carbons overall). A monounsaturated fatty acid must have exactly one carbon–carbon unsaturation 
(the unsaturation may be a double bond; any triple bonds disqualify the molecule) along its longest continuous 
(aliphatic) carbon chain.

This revised implementation:
  • Identifies the terminal carboxylic acid group using a SMARTS pattern.
  • Instead of rejecting molecules merely because some atoms occur in rings, it extracts the acyl chain by starting 
    from the sole carbon neighbour of the acid carbon and following only non‐ring C–C bonds.
  • Uses a DFS to find the longest continuous chain from that starting point.
  • Counts bond types along that path and accepts the molecule only if exactly one double bond (and zero triple bonds) is found.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES string.
      2. Finds a carboxylic acid group (using SMARTS, matching both protonated and deprotonated forms).
      3. Verifies that the acid carbon is terminal – that is, it has exactly one carbon neighbour.
      4. From that carbon neighbour (chain_start), constructs a graph of all carbon atoms connected by bonds 
         that are NOT in any ring (i.e. in the fatty acyl chain).
      5. Uses DFS to extract the longest continuous (acyclic) carbon chain starting from chain_start.
      6. Counts unsaturations (double bonds) and rejects any triple bonds along the chain.
      7. Finally, classifies the molecule as a monounsaturated fatty acid if the chain is sufficiently long 
         and has exactly one unsaturation.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        (bool, str): True and a success reason if molecule meets the criteria;
                     False and an explanation otherwise.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Locate a carboxylic acid group.
    # This pattern matches both protonated and deprotonated forms.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Molecule does not contain a carboxylic acid group, thus not a fatty acid"
    
    # Assume the first match is our target acid group.
    acid_match = acid_matches[0]
    acid_idx = acid_match[0]  # acid carbon is first atom in the pattern
    acid_atom = mol.GetAtomWithIdx(acid_idx)
    
    # The acid carbon should be terminal: it must have exactly one carbon neighbour.
    carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Acid group is not terminal; acid carbon does not have exactly one carbon neighbour"
    chain_start_atom = carbon_neighbors[0]
    
    # Ensure that the molecule overall has enough carbons to be a fatty acid.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 5:
        return False, "Not enough carbon atoms to be considered a fatty acid"
    
    # Build a graph for the continuous carbon chain.
    # Nodes are all carbon atoms (by their indices); edges exist only for carbon–carbon bonds that are NOT in any ring.
    # This allows us to traverse only the acyclic (aliphatic) portion of the molecule.
    graph = {}
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    for idx in carbon_indices:
        graph[idx] = []
    # Go through all bonds and only add connections for carbon atoms if the bond is not in a ring.
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if a.GetAtomicNum() == 6 and b.GetAtomicNum() == 6:
            if not bond.IsInRing():  # only use non-ring bonds for the fatty acyl chain
                graph[a.GetIdx()].append(b.GetIdx())
                graph[b.GetIdx()].append(a.GetIdx())
    
    # Use DFS to find the longest chain starting from chain_start_atom.
    longest_path = []
    def dfs(current, path, visited):
        nonlocal longest_path
        if len(path) > len(longest_path):
            longest_path = path.copy()
        for neighbor in graph[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                path.append(neighbor)
                dfs(neighbor, path, visited)
                path.pop()
                visited.remove(neighbor)
    
    start_idx = chain_start_atom.GetIdx()
    if start_idx not in graph:
        return False, "Chain-start atom is not present in the carbon graph"
    dfs(start_idx, [start_idx], {start_idx})
    
    # Require a minimum chain length (excluding the acid carbon).
    if len(longest_path) < 2:
        return False, "Fatty acid chain too short"
    
    # Count unsaturations along the consecutive bonds in the longest chain.
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
    
    # Final message: our chain (plus the terminal carboxyl group) qualifies as a monounsaturated fatty acid.
    return True, ("Molecule contains a terminal carboxylic acid group and an acyl carbon chain "
                  "with exactly one unsaturation (double bond) and no triple bonds, classifying it as a monounsaturated fatty acid")

# Example usage:
# if __name__ == '__main__':
#     test_smiles = "CCCCCCCC\\C=C/CCCCCCC(O)=O"  # cis-8-heptadecenoic acid
#     result, reason = is_monounsaturated_fatty_acid(test_smiles)
#     print(result, reason)