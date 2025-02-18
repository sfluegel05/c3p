"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA

A branched-chain fatty acyl-CoA is defined as a fatty acyl-CoA resulting 
from the condensation of the thiol group of coenzyme A with the carboxy group 
of a branched-chain fatty acid.

This updated program:
  1. Rejects molecules that appear deprotonated.
  2. Checks for a CoA moiety using an approximate SMARTS match.
  3. Looks for thioester groups (the acyl-CoA linkage).
  4. For each thioester candidate:
       a. Finds the carbonyl carbon and then identifies the acyl side (neighbors that are carbons, excluding the oxygen of the carbonyl).
       b. If more than one candidate acyl root is detected at the carbonyl, we immediately treat this as branching.
       c. Otherwise, starting at the acyl root, we collect all connected carbon atoms (via C–C bonds) while NOT going back into the carbonyl.
       d. We require that the candidate acyl fragment is acyclic and of minimum length.
       e. We then compute the “longest chain” (the diameter) within this induced subgraph. In a “linear” fragment the total number of atoms equals the longest chain length. Extra atoms (i.e. total > longest chain) indicate a branch.
  5. Returns True (with explanation) if any candidate acyl fragment shows extra carbons (branch) or if multiple acyl roots are found.

If none of the thioester candidates shows a branched acyl chain, the function returns False.
"""

from rdkit import Chem
from collections import deque

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        Tuple[bool, str]: A boolean indicating if the molecule is a branched-chain fatty acyl-CoA 
                          and a string with the reason for classification.
    """
    # 1. Reject molecules with explicit deprotonation markers.
    if "[O-]" in smiles:
        return False, "Molecule appears deprotonated ([O-] detected); expected neutral CoA form."
    
    # Parse molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 2. Check for a CoA moiety.
    # (This is an approximate SMARTS for part of the CoA fragment.)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found."
    
    # 3. Identify thioester groups. 
    # Here we look for a carbonyl carbon (C with =O) bonded to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found."
    
    # Helper: Given a starting carbon (acyl root) and excluding the carbonyl,
    # perform a DFS that collects all carbons connected via C–C bonds.
    def get_acyl_fragment(root_idx, carbonyl_idx):
        visited = set()
        stack = [root_idx]
        while stack:
            curr = stack.pop()
            if curr in visited:
                continue
            visited.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Do not traverse back into the carbonyl carbon
                if nbr_idx == carbonyl_idx:
                    continue
                if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                    stack.append(nbr_idx)
        return visited

    # Helper: Compute the longest chain length (diameter) in the induced subgraph.
    def longest_chain_length(fragment_indices):
        # Build a simple undirected graph for the fragment:
        graph = {idx: [] for idx in fragment_indices}
        for idx in fragment_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in fragment_indices:
                    if nbr_idx not in graph[idx]:
                        graph[idx].append(nbr_idx)
        # Breadth-first search helper.
        def bfs(start):
            distances = {start: 0}
            queue = deque([start])
            farthest = start
            while queue:
                current = queue.popleft()
                for nbr in graph[current]:
                    if nbr not in distances:
                        distances[nbr] = distances[current] + 1
                        queue.append(nbr)
                        if distances[nbr] > distances[farthest]:
                            farthest = nbr
            return farthest, distances[farthest]
        # Pick an arbitrary starting node.
        start = next(iter(fragment_indices))
        far_node, _ = bfs(start)
        _, diameter = bfs(far_node)
        return diameter + 1  # number of nodes in longest path
    
    # A minimal number of carbons in a valid fatty acyl fragment.
    MIN_FRAGMENT_LENGTH = 3

    # 4. For each thioester match, try to find and assess the acyl fragment.
    for match in thioester_matches:
        # In the SMARTS "[#6](=O)S", match[0] is the carbonyl carbon, match[1] is the =O,
        # and match[2] is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify candidate acyl roots from the carbonyl's neighbors (exclude oxygen and sulfur).
        candidate_acyl_roots = []
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in (match[1], sulfur_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                candidate_acyl_roots.append(nbr_idx)
        
        # If more than one candidate acyl root, this signals an immediate branch.
        if len(candidate_acyl_roots) > 1:
            explanation = ("Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety. "
                           "Multiple acyl carbon neighbors at the linkage indicate branching.")
            return True, explanation
        
        # If no candidate, move to next thioester candidate.
        if not candidate_acyl_roots:
            continue
        
        # With a single candidate acyl root, get the entire acyl fragment.
        acyl_root = candidate_acyl_roots[0]
        fragment = get_acyl_fragment(acyl_root, carbonyl_idx)
        
        # Require minimum fragment size.
        if len(fragment) < MIN_FRAGMENT_LENGTH:
            continue
        
        # Ensure the fragment does not contain ring atoms.
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in fragment):
            continue
        
        # Compute the longest linear chain (the diameter) in the fragment.
        longest_chain = longest_chain_length(fragment)
        # For a linear (unbranched) chain, total carbons == longest chain length.
        if len(fragment) > longest_chain:
            explanation = ("Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety. "
                           "Branching detected: total acyl fragment carbons (%d) exceed the longest chain length (%d)."
                           % (len(fragment), longest_chain))
            return True, explanation

    # If no thioester candidate shows a branched acyl chain:
    return False, "No branched acyl chain detected in any thioester candidate."


# Example usage (you can run this module as a script for a quick test):
if __name__ == "__main__":
    # Test with two examples:
    # 1. A branched one: 2-methylbutanoyl-CoA (should be True)
    smiles1 = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(smiles1)
    print("2-methylbutanoyl-CoA:", result, reason)
    
    # 2. A molecule that lacks a branched acyl chain (example: a linear fatty acyl-CoA analogue)
    smiles2 = "CCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(smiles2)
    print("Linear fatty acyl-CoA:", result, reason)