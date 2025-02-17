"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: Any branched-chain fatty acid containing methyl branches only.
This program defines the function is_methyl_branched_fatty_acid that takes a SMILES string as input,
and returns a bool along with a string explaining the classification.
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid.
    
    A methyl-branched fatty acid is defined as a fatty acid (i.e. containing a carboxylic acid group)
    with a long carbon chain that has one or more branches, 
    where each branch consists solely of a methyl (CH3) group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation supporting the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid functional group.
    # The SMARTS "C(=O)[O;H,-]" will match a carbonyl connected to an -OH (or its deprotonated form).
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Identify the carboxyl carbon (first atom in the matching substructure)
    matches = mol.GetSubstructMatches(acid_pattern)
    carboxyl_atom_idx = matches[0][0]
    
    # Build a carbon-only subgraph.
    # Consider only atoms with atomic number 6 and bonds between them.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_idx_set = set(atom.GetIdx() for atom in carbon_atoms)
    carbon_adj = {idx: [] for idx in carbon_idx_set}
    for idx in carbon_idx_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                n_idx = nbr.GetIdx()
                # Add neighbor if it is also a carbon.
                if n_idx in carbon_idx_set:
                    carbon_adj[idx].append(n_idx)
    
    # Use a depth-first search starting at the carboxyl carbon to determine the longest carbon chain.
    longest_chain = []
    def dfs(current, path, visited):
        nonlocal longest_chain
        # Update the longest chain found so far.
        if len(path) > len(longest_chain):
            longest_chain = path.copy()
        for neighbor in carbon_adj[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                path.append(neighbor)
                dfs(neighbor, path, visited)
                path.pop()
                visited.remove(neighbor)
    
    # Ensure the carboxyl carbon is in our carbon set (it should be, as it is carbon).
    if carboxyl_atom_idx not in carbon_idx_set:
        return False, "Carboxyl carbon not found in carbon backbone"
    dfs(carboxyl_atom_idx, [carboxyl_atom_idx], {carboxyl_atom_idx})
    main_chain_set = set(longest_chain)
    
    # Check for branching: each branch is a carbon attached to a main chain carbon that is not part of the main chain.
    # We require that each such branch is a methyl group (i.e. it has only one carbon neighbor and exactly 3 hydrogens).
    branch_count = 0
    for main_idx in main_chain_set:
        main_atom = mol.GetAtomWithIdx(main_idx)
        for nbr in main_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in main_chain_set:
                # For a methyl branch, the sole carbon neighbor should be the main chain atom.
                carbon_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 6]
                if len(carbon_neighbors) != 1:
                    return False, f"Found a branch attached to atom index {main_idx} that has more than one carbon neighbor (i.e. not a simple methyl group)"
                # Check the hydrogen count. (Implicit hydrogens are taken into account with GetTotalNumHs.)
                if nbr.GetTotalNumHs() != 3:
                    return False, f"Branch attached to atom index {main_idx} is not a methyl (expected 3 hydrogens, got {nbr.GetTotalNumHs()})"
                branch_count += 1

    if branch_count == 0:
        return False, "No branches found; the fatty acid appears to be a straight-chain acid"
    
    # (Optional) Verify that the main chain is long enough to be considered a fatty acid.
    if len(main_chain_set) < 5:
        return False, "Main carbon chain is too short to be a fatty acid"
    
    return True, f"Methyl-branched fatty acid with {branch_count} methyl branch(es) on a main chain of length {len(main_chain_set)}"
    
# Example usage (if run as a script, one could test on an example SMILES):
if __name__ == '__main__':
    test_smiles = "CC(C)CCCCCCCCCCCCCCCCC(O)=O"  # 18-methylnonadecanoic acid
    result, reason = is_methyl_branched_fatty_acid(test_smiles)
    print(result, reason)