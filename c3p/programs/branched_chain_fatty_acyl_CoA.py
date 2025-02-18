"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA

A branched-chain fatty acyl-CoA is defined as a fatty acyl-CoA in which the acyl (fatty acid) moiety
contains at least one branch (i.e. a substituent off the “main” carbon chain). In practice, we must detect
(i) a thioester linkage –C(=O)S– joining the acyl fragment to the Coenzyme A moiety, 
(ii) the presence of a CoA substructure, and 
(iii) branching in the fatty acyl fragment.
To determine branching we:
    • first extract a candidate acyl fragment starting from the carbon attached to the carbonyl.
    • insist that the fragment is acyclic (no ring atoms).
    • use two complementary tests: (a) check if the acyl root shows extra carbon neighbors,
      and (b) compare the number of carbons in the acyl fragment to the length of its longest acyclic carbon chain.
If either test indicates extra branches, we classify the molecule as a branched-chain fatty acyl-CoA.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

    Steps:
      1. Parse the SMILES string.
      2. Verify that the molecule contains a CoA motif (using an approximate SMARTS).
      3. Identify a thioester group with the SMARTS pattern "[#6](=O)S".
      4. For each thioester occurrence:
            a. Identify the acyl (fatty acid) fragment. We take the neighbor (carbon) of the 
               carbonyl (excluding the oxygen and sulfur) as the “acyl root” and traverse 
               through carbon–carbon bonds.
            b. Reject the candidate if any of its atoms is in a ring (we expect a linear chain).
            c. Test for branching via two complementary ideas:
                   (i) Check if at the acyl root extra carbon neighbors occur beyond the expected one (or one plus the carbonyl).
                  (ii) Compute the longest simple (acyclic) carbon chain in the fragment. In a linear chain,
                       this length equals the total number of carbons; if shorter then branching exists.
      5. Return True with an appropriate message if any thioester yields a branched acyl fragment.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if the molecule is classified as branched-chain fatty acyl-CoA, False otherwise.
       str: Reason for the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for the CoA moiety.
    # This SMARTS roughly matches a core substructure found in most CoA esters.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"
    
    # Step 2: Identify a thioester group.
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found"
    
    # Helper function: starting from the acyl root, traverse the acyl fragment.
    # Do not go back to the carbonyl (donor) atom.
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
                # Do not return to the carbonyl atom
                if nbr_idx == carbonyl_idx:
                    continue
                # Only follow carbon neighbors
                if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                    frontier.append(nbr_idx)
        return visited

    # Helper function: compute the longest simple (acyclic) path length (number of atoms)
    # within the set of atom indices provided. (Uses DFS on the subgraph defined by the acyl fragment.)
    def longest_path_length(atom_indices):
        sub_atoms = list(atom_indices)
        n = len(sub_atoms)
        # Map original idx to local index
        local_idx = {a: i for i, a in enumerate(sub_atoms)}
        # Create adjacency list for atoms in the fragment
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
        best = 0
        def dfs(u, visited, length):
            nonlocal best
            best = max(best, length)
            for v in adj[u]:
                if v not in visited:
                    dfs(v, visited | {v}, length+1)
        for i in range(n):
            dfs(i, {i}, 1)
        return best

    # Now process each thioester match. We require that at least one thioester yields a branched acyl fragment.
    branch_detected = False
    explanation = ""
    for match in thioester_matches:
        # In the SMARTS "[#6](=O)S", match[0]=carbonyl carbon, match[1]=oxygen, match[2]=sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Find the neighbor from carbonyl that is not oxygen or sulfur. This neighbor is the acyl root.
        acyl_root = None
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in (match[1], sulfur_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_root = nbr_idx
                break
        if acyl_root is None:
            continue  # try next thioester match if available

        # Extract the acyl fragment.
        acyl_fragment = get_acyl_fragment(acyl_root, carbonyl_idx)
        if len(acyl_fragment) < 2:
            continue  # not enough carbons in the fatty acyl chain

        # NEW: Reject if any atom in the acyl fragment is part of a ring.
        if any(mol.GetAtomWithIdx(a).IsInRing() for a in acyl_fragment):
            # This candidate is likely not a simple fatty acyl chain.
            continue

        # Test (a): Check acyl_root for branching.
        # In a simple unbranched chain the acyl root should have one carbon neighbor (inside the chain)
        # aside from the connection to the carbonyl.
        acyl_root_atom = mol.GetAtomWithIdx(acyl_root)
        carbon_count = 0
        for nbr in acyl_root_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx == carbonyl_idx or nbr_idx in acyl_fragment:
                carbon_count += 1
        # For a linear chain this count should be at most 2 (one from the carbonyl and one acyl neighbor).
        if carbon_count > 2:
            branch_detected = True
            explanation = "Contains a branched acyl chain (extra carbon neighbor at acyl root)."
            break

        # Test (b): Compare number of carbons in the acyl fragment to longest simple chain.
        lp = longest_path_length(acyl_fragment)
        if len(acyl_fragment) > lp:
            branch_detected = True
            explanation = "Contains a branched acyl chain (total carbons exceed longest chain length)."
            break

    if not branch_detected:
        return False, "No branched fatty acyl chain detected"
    
    return True, f"Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety. {explanation}"


# Example usage:
if __name__ == "__main__":
    # Test with a known branched-chain fatty acyl-CoA example: 2-methylbutanoyl-CoA
    example_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(example_smiles)
    print(result, reason)