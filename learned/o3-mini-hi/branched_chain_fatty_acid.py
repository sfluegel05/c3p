"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid
A branched-chain fatty acid (BCFA) is defined as any fatty acid with one or more alkyl substituents
on its parent hydrocarbon chain (usually a carboxylic acid with a long alkyl chain that contains at least one branch).
"""

from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid (BCFA) based on its SMILES string.
    
    Approach:
      1. Check for a carboxylic acid group in the molecule using a SMARTS pattern.
         (A fatty acid must contain a -C(=O)OH group.)
      2. Build a graph of all carbon atoms (atomic number 6) and their carbon–carbon connectivity.
      3. From the carboxyl carbon (the carbon of the -COOH group), determine the longest continuous chain
         (assumed to be the fatty acyl chain) via a depth-first search.
      4. For each atom in this main chain, check if there is any carbon neighbor that is not on the main chain;
         if yes, then an alkyl substituent (branch) is present.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule qualifies as a branched-chain fatty acid, False otherwise.
        str: A reason that explains the decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check for a carboxylic acid group.
    # The SMARTS pattern [CX3](=O)[OH] matches a typical carboxyl group.
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OH]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid."
    
    # For simplicity, assume the first match corresponds to the fatty acid’s carboxyl group.
    # The carboxyl carbon is the first atom in the match.
    carboxyl_idx = acid_matches[0][0]
    
    # Build a subgraph of carbon atoms.
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 3:
        return False, "Too few carbon atoms to be a fatty acid."
    
    # Build connectivity (neighbors) among carbon atoms.
    carbon_graph = {idx: [] for idx in carbon_atoms}
    for idx in carbon_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_graph[idx].append(neighbor.GetIdx())
    
    # Define a DFS function to find the longest simple path (in terms of atom count)
    def dfs(current, visited):
        best_path = [current]
        for neighbor in carbon_graph[current]:
            if neighbor not in visited:
                path = [current] + dfs(neighbor, visited | {neighbor})
                if len(path) > len(best_path):
                    best_path = path
        return best_path
    
    # Ensure the carboxyl carbon is in our carbon_atoms (it should be)
    if carboxyl_idx not in carbon_graph:
        return False, "Carboxyl carbon not found in carbon subgraph."
    
    # Find the longest chain (main chain) starting from the carboxyl carbon.
    main_chain = dfs(carboxyl_idx, {carboxyl_idx})
    if len(main_chain) < 3:
        return False, "The fatty acyl chain appears too short."
    
    # Now, for each atom in the main chain, check if there is a carbon substituent that is not part of the main chain.
    branch_found = False
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Look at carbon neighbors (in the whole molecule)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                neigh_idx = neighbor.GetIdx()
                if neigh_idx not in main_chain:
                    branch_found = True
                    break
        if branch_found:
            break
    
    if not branch_found:
        return False, "No alkyl substituents (branch) found on the fatty acyl chain."
    
    # Optionally, we can add a check that the molecule is mostly acyclic in the fatty acid part.
    # For many fatty acids the chain is linear aside from the branch.
    # (Additional criteria may be included based on chain length or saturation if desired.)
    
    return True, "Contains a carboxylic acid group with a fatty acyl chain that has an alkyl branch."

# Example usage:
if __name__ == "__main__":
    # Test with an example: isobutyric acid (CC(C)C(O)=O) should be a branched-chain fatty acid.
    test_smiles = "CC(C)C(O)=O"  
    result, reason = is_branched_chain_fatty_acid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")