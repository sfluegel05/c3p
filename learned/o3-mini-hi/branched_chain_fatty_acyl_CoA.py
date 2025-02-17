"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA

A branched-chain fatty acyl-CoA is defined as a fatty acyl-CoA in which the acyl (fatty acid) moiety
contains at least one branch (i.e. a substituent off the “main” carbon chain). In practice, we must detect
(i) a thioester linkage –C(=O)S– joining the acyl fragment to the Coenzyme A moiety, (ii) the presence of
a CoA substructure, and (iii) branching in the fatty acyl fragment. Because the acyl fragment is all the carbon
atoms connected to the carbonyl (other than the carbonyl itself), a linear chain will have a longest continuous path
equal in length to its total number of carbons. In a branched chain that number will be higher.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

    The algorithm performs the following steps:
      1. Parse the SMILES and check for a Coenzyme A (CoA) fragment – we look for a characteristic pattern.
      2. Identify a thioester group using the SMARTS pattern "[#6](=O)S", i.e. a carbonyl carbon attached to sulfur.
      3. For each thioester found, determine the acyl (fatty acid) fragment. This is done by taking the carbon
         neighbor of the carbonyl (other than the oxygen and the sulfur) as the “acyl root” and then traversing
         via carbon–carbon bonds (avoiding going back into the carbonyl).
      4. To determine whether the acyl fragment is branched, we use two complementary ideas:
             (a) For every carbon (in the acyl fragment) that is directly bonded to more than two other carbons,
                 count such degree after adding the carbonyl bond to the acyl root. In a linear chain, an internal
                 carbon will have only two carbon neighbors.
             (b) Compute the longest simple (acyclic) carbon path within the acyl fragment and compare its length
                 to the total number of carbons in the fragment. For a linear chain these numbers are equal;
                 for a branched chain the total count exceeds the longest path.
         If either of these tests indicates branching, the molecule is classified as a branched-chain fatty acyl-CoA.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a branched-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for the CoA moiety.
    # Use a SMARTS pattern that matches a core substructure present in most CoA esters.
    # (Note: CoA is large and variable; this pattern is an approximation.)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"
    
    # Step 2: Identify the thioester group (the acyl-CoA linkage).
    # We search for a carbon with a double-bonded O and a single bond to S.
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found"
    
    # Helper: traverse the acyl chain starting from a given atom (acyl root)
    # Exclude going back to the carbonyl carbon.
    def get_acyl_fragment(root_idx, carbonyl_idx):
        """Return the set of carbon atom indices in the acyl fragment (excluding the carbonyl carbon)."""
        visited = set()
        frontier = [root_idx]
        while frontier:
            current = frontier.pop()
            if current in visited:
                continue
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Do not return to the carbonyl carbon
                if nbr_idx == carbonyl_idx:
                    continue
                # Only follow carbon neighbors
                if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                    frontier.append(nbr_idx)
        return visited

    # Helper: compute the longest simple (acyclic) path length (number of atoms) within the subgraph defined by indices.
    def longest_path_length(atom_indices):
        # Build a subgraph (list of atoms) and record bonds among them
        sub_atoms = list(atom_indices)
        n = len(sub_atoms)
        # Map original idx to local index
        local_idx = {a: i for i, a in enumerate(sub_atoms)}
        # Build an adjacency list for the subgraph.
        adj = [[] for _ in range(n)]
        for a in sub_atoms:
            atom = mol.GetAtomWithIdx(a)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in atom_indices:
                    i = local_idx[a]
                    j = local_idx[nbr_idx]
                    if j not in adj[i]:
                        adj[i].append(j)
                    if i not in adj[j]:
                        adj[j].append(i)
        # For each vertex, run a DFS to compute the longest simple path.
        best = 0
        def dfs(u, prev, length, visited):
            nonlocal best
            best = max(best, length)
            for v in adj[u]:
                if v == prev or v in visited:
                    continue
                dfs(v, u, length+1, visited | {v})
        for i in range(n):
            dfs(i, -1, 1, {i})
        return best

    # Now process each thioester match.
    branch_detected = False
    for match in thioester_matches:
        # In the SMARTS "[#6](=O)S", match[0]=carbonyl carbon, match[1]=oxygen, match[2]=sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Find the neighbor to the carbonyl that is not the oxygen or the sulfur; this is the acyl root.
        acyl_root = None
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in (match[1], sulfur_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_root = nbr_idx
                break
        if acyl_root is None:
            continue  # try next thioester match if any
        
        # Get the set of carbon atoms in the acyl fragment.
        acyl_fragment = get_acyl_fragment(acyl_root, carbonyl_idx)
        if len(acyl_fragment) < 2:
            continue  # not enough carbons to be a fatty acyl chain
        
        # Test (a): Check if the acyl root (which is attached to the carbonyl) shows branching.
        # For acyl_root, count carbon neighbors that are either in the acyl fragment or (if attached) the carbonyl.
        acyl_root_atom = mol.GetAtomWithIdx(acyl_root)
        count_neighbors = 0
        for nbr in acyl_root_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            # Include the carbonyl if attached.
            if nbr_idx == carbonyl_idx or nbr_idx in acyl_fragment:
                count_neighbors += 1
        # In an unbranched (primary) chain the acyl root will have only 1 carbon neighbor from the acyl fragment (or 1 + carbonyl).
        if count_neighbors > 2:
            branch_detected = True
            break
        
        # Test (b): Compare the total number of carbons in the acyl fragment with the length of its longest path.
        # In a linear chain these numbers are equal; a branched chain has extra side groups.
        lp = longest_path_length(acyl_fragment)
        if len(acyl_fragment) > lp:
            branch_detected = True
            break

    if not branch_detected:
        return False, "No branched fatty acyl chain detected"
    
    return True, "Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety"


# Example usage:
if __name__ == "__main__":
    # For example, test with 2-methylbutanoyl-CoA (a known branched-chain fatty acyl-CoA)
    example_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(example_smiles)
    print(result, reason)