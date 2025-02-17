"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: Any branched-chain fatty acid containing methyl branches only.
This program defines the function is_methyl_branched_fatty_acid that takes a SMILES string as input,
and returns a bool along with a reason string for the classification.
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid.

    A methyl-branched fatty acid is defined as a fatty acid containing a carboxylic acid group 
    with a long carbon chain that has one or more substituents, 
    where every branch is a simple methyl (-CH3) group.
    
    This function first confirms that a carboxylic acid group is present.
    Then it builds a subgraph of carbons and uses DFS to find the longest carbon chain starting
    from the carboxyl carbon. Once that main chain is established, it checks that any extra carbon
    (i.e. any branch off the main chain) is a methyl group (attached only to the main chain with exactly 3 hydrogens).
    Finally, if there are extra carbons in the molecule beyond the main chain, they must be exactly those
    methyl branches.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation supporting the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid: pattern matches carbonyl bonded to hydroxyl (or its deprotonated form)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Use the first match; assume the first atom is the carboxyl carbon.
    carboxyl_idx = acid_matches[0][0]
    
    # Build a carbon-only graph: a dictionary mapping carbon atom index to a list of neighboring carbon indices.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found in molecule"
    carbon_idx_set = set(atom.GetIdx() for atom in carbon_atoms)
    carbon_adj = {idx: [] for idx in carbon_idx_set}
    for idx in carbon_idx_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                n_idx = nbr.GetIdx()
                if n_idx in carbon_idx_set:
                    carbon_adj[idx].append(n_idx)
                    
    # Check that the carboxyl carbon is in our carbon set.
    if carboxyl_idx not in carbon_idx_set:
        return False, "Carboxyl carbon not found in carbon backbone"
    
    # Find the longest simple carbon chain starting from the carboxyl carbon.
    longest_chain = []
    def dfs(current, path, visited):
        nonlocal longest_chain
        # Update the longest chain if current path is longer.
        if len(path) > len(longest_chain):
            longest_chain = path.copy()
        for neighbor in carbon_adj.get(current, []):
            if neighbor not in visited:
                visited.add(neighbor)
                path.append(neighbor)
                dfs(neighbor, path, visited)
                path.pop()
                visited.remove(neighbor)
                
    dfs(carboxyl_idx, [carboxyl_idx], {carboxyl_idx})
    main_chain = longest_chain
    main_chain_set = set(main_chain)
    
    # Enforce a minimum main chain length to qualify as a fatty acid.
    if len(main_chain) < 5:
        return False, "Main carbon chain is too short to be a fatty acid"
    
    # Helper: Determine whether an atom is a methyl group attached to the main chain.
    def is_simple_methyl(atom, attached_to_idx):
        # It must be a carbon, not in the main chain.
        if atom.GetAtomicNum() != 6:
            return False
        # Among its neighbors, if we only see the one it is attached to (i.e. degree one in the carbon graph) then it is methyl.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1 or carbon_neighbors[0].GetIdx() != attached_to_idx:
            return False
        # Check that it has exactly 3 hydrogens (implicit+explicit).
        if atom.GetTotalNumHs() != 3:
            return False
        return True
    
    branch_count = 0
    # For each carbon in the main chain, check neighbors not in the main chain.
    for main_idx in main_chain:
        main_atom = mol.GetAtomWithIdx(main_idx)
        for nbr in main_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in main_chain_set:
                # It must be a simple methyl.
                if not is_simple_methyl(nbr, main_idx):
                    return False, f"Found a branch attached to atom index {main_idx} that is not a simple methyl group"
                branch_count += 1

    if branch_count == 0:
        return False, "No branches found; the fatty acid appears to be a straight-chain acid"
    
    # Count all carbons in molecule and subtract those in main chain.
    total_carbons = len(carbon_idx_set)
    extra_carbons = total_carbons - len(main_chain_set)
    if extra_carbons != branch_count:
        return False, ("Extra carbon atoms beyond the methyl branches were found "
                       f"({extra_carbons} extra vs {branch_count} methyl branch(es)).")
    
    # Build classification message.
    return True, (f"CORRECT Methyl-branched fatty acid with {branch_count} methyl branch(es) "
                  f"on a main chain of length {len(main_chain_set)}")

# Example usage (if run as a script, one could test on an example SMILES):
if __name__ == '__main__':
    # Test on one of the given examples: 18-methylnonadecanoic acid
    test_smiles = "CC(C)CCCCCCCCCCCCCCCCC(O)=O"
    result, reason = is_methyl_branched_fatty_acid(test_smiles)
    print(result, reason)