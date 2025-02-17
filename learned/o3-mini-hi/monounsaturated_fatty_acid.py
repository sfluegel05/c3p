"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated Fatty Acid

Definition: A fatty acid is defined to have a terminal carboxylic acid group and an acyclic carbon chain 
(with at least 5 carbons overall). A monounsaturated fatty acid must have exactly one carbon–carbon unsaturation 
(the unsaturation may be a double bond; any triple bonds disqualify the molecule) along its longest continuous 
(aliphatic) carbon chain.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    
    This implementation:
      1. Parses the SMILES string.
      2. Rejects cyclic molecules.
      3. Identifies a terminal carboxylic acid group via a SMARTS pattern that catches both protonated 
         and deprotonated forms.
      4. Checks that the acid carbon (of the carboxyl group) is terminal (i.e. it has exactly one carbon neighbor).
      5. Builds a carbon-only graph (excluding the acid carbon) and 
         extracts the longest continuous carbon chain starting from the neighbor of the acid group.
      6. Counts the unsaturations (bonds whose bond type is DOUBLE) along that chain and immediately rejects any 
         triple bonds.
      7. Finally, classifies the molecule as a monounsaturated fatty acid only if exactly one double 
         bond is found in the chain.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        (bool, str): Tuple with True and a reason message if the classification is positive, 
                     else False and an explanation.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject cyclic molecules.
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule is cyclic, which is not typical for a fatty acid chain"
    
    # Identify a carboxylic acid group (matches both protonated and deprotonated forms)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Molecule does not contain a carboxylic acid group, thus not a fatty acid"
    
    # Take the first acid match. We assume the acid carbon is the first atom in the match.
    acid_match = acid_matches[0]
    acid_idx = acid_match[0]
    acid_atom = mol.GetAtomWithIdx(acid_idx)
    
    # The carboxylic acid carbon should be terminal: it must have exactly one neighbor that is a carbon.
    carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Acid group is not terminal; the acid carbon has an unexpected bonding pattern"
    
    # This neighbor initiates our fatty acid chain.
    chain_start = carbon_neighbors[0].GetIdx()
    
    # Ensure the molecule has a sufficient total number of carbons (including the acid carbon) to be a fatty acid.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 5:
        return False, "Molecule does not have enough carbon atoms to be considered a fatty acid"
    
    # Build a graph of all carbon atoms (as nodes) ignoring the acid carbon.
    # We only consider bonds between carbons.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() 
                        if atom.GetAtomicNum() == 6 and atom.GetIdx() != acid_idx]
    graph = {i: [] for i in carbon_indices}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Only count bonds between carbons that are in our graph.
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and 
            a1.GetIdx() in graph and a2.GetIdx() in graph):
            graph[a1.GetIdx()].append(a2.GetIdx())
            graph[a2.GetIdx()].append(a1.GetIdx())
    
    # Use DFS to find the longest simple path (greatest number of atoms) starting from chain_start.
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
                
    if chain_start not in graph:
        return False, "Chain-start atom is not present in the carbon chain graph"
    dfs(chain_start, [chain_start], {chain_start})
    
    # Count unsaturations in the extracted longest chain.
    # We count only bonds along the consecutive atoms of the path.
    double_count = 0
    triple_count = 0
    for i in range(len(longest_path) - 1):
        a_idx = longest_path[i]
        b_idx = longest_path[i+1]
        bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
        if bond is None:
            continue
        btype = bond.GetBondType()
        if btype == Chem.BondType.DOUBLE:
            double_count += 1
        elif btype == Chem.BondType.TRIPLE:
            triple_count += 1
    
    # Immediately reject molecules with any triple bonds in the chain.
    if triple_count > 0:
        return False, "Molecule has a triple bond in the fatty acid chain"
    
    # The classification requires exactly one double (unsaturation) bond.
    if double_count != 1:
        return False, f"Fatty acid chain unsaturation count is {double_count} rather than the required 1"
    
    # Optionally check that the chain is of a reasonable length (here at least 3 carbons in the chain beyond the acid group).
    if len(longest_path) < 3:
        return False, "Carbon chain too short to be a fatty acid"
    
    return True, ("Molecule contains a terminal carboxylic acid group and an acyclic carbon chain with exactly one "
                  "unsaturation (double bond) and no triple bonds, classifying it as a monounsaturated fatty acid")

# Example usage:
# if __name__ == "__main__":
#     # Test with cis-8-heptadecenoic acid:
#     test_smiles = "CCCCCCCC\\C=C/CCCCCCC(O)=O"
#     result, reason = is_monounsaturated_fatty_acid(test_smiles)
#     print(result, reason)