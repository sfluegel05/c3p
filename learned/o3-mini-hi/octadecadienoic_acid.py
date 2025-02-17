"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid – any straight-chain C18 polyunsaturated fatty acid
with exactly 2 non-aromatic C=C double bonds and a terminal carboxylic acid group.
"""

from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string)
    is an octadecadienoic acid. The definition enforced here is:
      • the molecule possesses a carboxylic acid group at one terminus,
      • the molecule’s longest carbon chain (i.e. the main chain) has 18 carbon atoms,
      • along that main chain there are exactly 2 non-aromatic C=C double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple. True with a positive reason if it is an octadecadienoic acid,
                     otherwise False with an explanatory reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid functionality using SMARTS.
    # We expect a terminal carboxylic acid: the acid carbon inside "[CX3](=O)[OX2H1]".
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing carboxylic acid functionality"

    # Build an undirected graph on all carbon atoms.
    # Nodes are atom indices for carbons and edges exist between two carbons that are chemically bonded.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found"
    carbon_graph = {}
    for atom in carbon_atoms:
        idx = atom.GetIdx()
        # Only add neighbors that are carbons.
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        carbon_graph[idx] = neighbors

    # Helper function to find the longest simple path (as a list of carbon indices)
    # from a given starting node in the carbon graph.
    def dfs(current, visited):
        best_path = [current]
        for nbr in carbon_graph.get(current, []):
            if nbr in visited:
                continue
            candidate = [current] + dfs(nbr, visited | {nbr})
            if len(candidate) > len(best_path):
                best_path = candidate
        return best_path

    # Find the overall longest carbon chain.
    longest_chain = []
    # Try starting a DFS from each carbon atom.
    for start in carbon_graph.keys():
        path = dfs(start, {start})
        if len(path) > len(longest_chain):
            longest_chain = path

    # For a straight-chain C18 fatty acid the main carbon chain must have exactly 18 carbons.
    if len(longest_chain) != 18:
        return False, f"Longest carbon chain has {len(longest_chain)} carbons instead of 18"

    # Check that the chain is linear.
    # In a simple path the connectivity within the chain is linear by construction,
    # but we also check that the path is not a detour from a branched network.
    # (We do not reject branches in the overall molecule; only the main chain must be unbranched.)
    # For each carbon in the longest chain, count how many neighbors are also in the chain.
    chain_set = set(longest_chain)
    for i, idx in enumerate(longest_chain):
        # Count neighbors in the chain
        in_chain_neighbors = sum(1 for nbr in carbon_graph[idx] if nbr in chain_set)
        # Endpoints should have exactly one neighbor in the chain,
        # internal carbons exactly 2.
        if i == 0 or i == len(longest_chain)-1:
            if in_chain_neighbors != 1:
                return False, f"Terminal carbon {idx} in main chain does not have exactly one chain neighbor"
        else:
            if in_chain_neighbors != 2:
                return False, f"Internal carbon {idx} in main chain does not have exactly two chain neighbors"
    
    # Check that the terminal carbon at one end is part of the carboxylic acid group.
    # The carboxyl carbon is the first atom in the SMARTS match.
    acid_carbons = {match[0] for match in acid_matches}
    if longest_chain[0] not in acid_carbons and longest_chain[-1] not in acid_carbons:
        return False, "Neither end of the main chain shows carboxylic acid functionality"
    
    # Count the number of non-aromatic double bonds along the main chain.
    double_bond_count = 0
    for i in range(len(longest_chain) - 1):
        bond = mol.GetBondBetweenAtoms(longest_chain[i], longest_chain[i+1])
        if bond is None:
            return False, f"Missing bond between carbons {longest_chain[i]} and {longest_chain[i+1]}"
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
            double_bond_count += 1
    
    if double_bond_count != 2:
        return False, f"Found {double_bond_count} non-aromatic C=C bonds along the main chain; exactly 2 are required"
    
    return True, ("Molecule is a straight-chain C18 fatty acid with exactly 2 non-aromatic C=C bonds "
                  "and a terminal carboxylic acid group")

# For testing (this block will run if the module is executed directly)
if __name__ == "__main__":
    # Sample test cases (including some examples from the provided dataset)
    test_cases = [
        # Expected true positives:
        ("CCCCCC\\C=C/C=C/[C@H](O)CCCCCCCC(O)=O", "9(R)-HODE"),
        ("OC(=O)CCCCCCC\\C=C/C=C\\CCCCCC", "9Z,11Z-octadecadienoic acid"),
        ("CC/C=C\\C/C=C\\CC(C(CCCCCCCC(O)=O)O)O", "9,10-DiHODE"),
        # Expected false positives or negatives:
        ("C(=C\\C/C=C\\CCCCCO)CCCCCCCC(=O)O", "18-hydroxylinoleic acid (should be rejected)"),
        ("OC(=O)CCCCC/C=C/C=C\\CCCCCCCC", "7-trans,9-cis-octadecadienoic acid (should be rejected)"),
        # A false negative example (extra carbon in a side chain)
        ("O(C(CCCCCCC(O)=O)/C=C/C(=O)CCCCCCCCC(O)=O", "(11E)-13-hydroxy-10-oxo-11-octadecenoic acid (should be rejected)"),
    ]
    
    for sm, name in test_cases:
        res, reason = is_octadecadienoic_acid(sm)
        status = "CORRECT" if res else "REJECTED"
        print(f"SMILES: {sm}\nNAME: {name}\nResult: {status}\nReason: {reason}\n")