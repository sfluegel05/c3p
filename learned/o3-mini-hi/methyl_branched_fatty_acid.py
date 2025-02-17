"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: Any branched‐chain fatty acid containing methyl branches only.
This function, is_methyl_branched_fatty_acid, takes a SMILES string as input and returns a boolean and a reason.
The approach:
  1. Parse molecule and check for a carboxylic acid group.
  2. Reject molecules with rings (non‐acyclic structures).
  3. Build a carbon-only graph.
  4. Enumerate all simple paths (starting at the acid carbon) as possible backbones.
  5. For each candidate backbone (which must be at least 5 carbons long) check that for every carbon in the chain the branch 
     (i.e. any adjacent carbons not in the backbone) is “methyl‐only”. For a branch we allow either a single CH3 group
     or a star‐shaped (one–level only) branch (as in isopropyl or tert‐butyl) where the attached “root” carbon’s neighbors 
     (all in the branch) are terminal methyl groups.
  6. If one candidate backbone is found, return True and an explanation that reports the number of branches (each branch counts as one)
     and the backbone length.
  
If no candidate is found, a False is returned along with a reason.
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Reject molecules with rings (fatty acids are expected acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a straight-chain fatty acid"
        
    # Look for a carboxylic acid group (we allow for deprotonated forms too)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    # choose the first carboxyl carbon (the carbonyl carbon)
    carboxyl_idx = acid_matches[0][0]
    
    # Build a carbon-only graph: dictionary mapping carbon atom index to set of neighboring carbon indices.
    carbon_idx_set = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum()==6}
    if not carbon_idx_set:
        return False, "No carbon atoms present"
    carbon_adj = {idx: set() for idx in carbon_idx_set}
    for idx in carbon_idx_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum()==6 and nbr.GetIdx() in carbon_idx_set:
                carbon_adj[idx].add(nbr.GetIdx())
                
    # Ensure that the carboxyl carbon is in the carbon graph.
    if carboxyl_idx not in carbon_idx_set:
        return False, "Carboxyl carbon not in carbon backbone"
        
    # Helper function: Check if a branch attached to a backbone atom is "methyl-only".
    # The branch is defined as the connected subgraph (in carbon graph) starting at branch_root,
    # staying completely outside the current backbone.
    # We allow either a single carbon (which must be CH3) or a one-level star (the root with only terminal CH3 groups).
    def valid_branch(branch_root, backbone_set):
        # Perform a BFS on the branch subgraph (only nodes not in backbone_set)
        from collections import deque
        q = deque()
        q.append((branch_root, 0))  # (node, distance from branch_root)
        branch_nodes = {branch_root: 0}
        while q:
            current, dist = q.popleft()
            # We do not allow depth > 1 (except for the root at depth 0).
            if dist > 1:
                return False, f"Branch starting at atom {branch_root} extends beyond allowed depth"
            # For each neighbor in the carbon graph that is not in the backbone:
            for nbr in carbon_adj[current]:
                if nbr in backbone_set:
                    continue
                if nbr not in branch_nodes:
                    branch_nodes[nbr] = dist+1
                    if dist+1 > 1:
                        return False, f"Branch starting at atom {branch_root} has chain length >1"
                    q.append((nbr, dist+1))
        # Now, if branch is a single carbon (only branch_root found)
        branch_atom = mol.GetAtomWithIdx(branch_root)
        if len(branch_nodes)==1:
            # It must be a methyl: exactly 3 hydrogens.
            if branch_atom.GetTotalNumHs() != 3:
                return False, f"Atom {branch_root} is not a CH3 group"
            return True, None
        else:
            # Multi-carbon branch: By construction, only two layers exist: the root (depth 0) and its immediate neighbors (depth 1).
            # Let d1 be all branch nodes with distance 1.
            d1 = [n for n,d in branch_nodes.items() if d==1]
            # Check that the root's hydrogen count equals 3 - (number of neighbors in branch)
            expected_root_H = 3 - len(d1)
            if branch_atom.GetTotalNumHs() != expected_root_H:
                return False, f"Branch root (atom {branch_root}) hydrogen count ({branch_atom.GetTotalNumHs()}) does not equal expected ({expected_root_H})"
            # Now, each neighbor in d1 must be a terminal CH3 (no further branch in branch graph, so no other carbon neighbor in branch)
            for nbr in d1:
                nbr_atom = mol.GetAtomWithIdx(nbr)
                # count how many neighbors are within branch_nodes (should be 1, linking back to branch root)
                cnt = 0
                for nn in carbon_adj[nbr]:
                    if nn in branch_nodes:
                        cnt += 1
                # Terminal methyl should be attached only to the root.
                if cnt != 1 or nbr_atom.GetTotalNumHs() != 3:
                    return False, f"Atom {nbr} in branch is not a terminal CH3 group"
            return True, None

    # Helper: Enumerate all simple paths (backbones) starting from the carboxyl carbon.
    # Because the carbon graph is acyclic (we rejected rings) we can use recursive DFS.
    all_paths = []
    def dfs_path(current, path):
        # Extend DFS by exploring all neighbors not already in path.
        extended = False
        for nbr in carbon_adj[current]:
            if nbr in path:
                continue
            dfs_path(nbr, path + [nbr])
            extended = True
        # If no extension, record the current path as a candidate leaf path.
        if not extended:
            all_paths.append(path)
    dfs_path(carboxyl_idx, [carboxyl_idx])
    
    # Among the candidate paths, choose the one with maximal length that also meets the branched criteria.
    valid_candidates = []
    for chain in all_paths:
        # Enforce a minimum chain length
        if len(chain) < 5:
            continue
        backbone_set = set(chain)
        branch_errors = []
        branch_count = 0
        # For every atom in the chain, check neighbors not in the backbone.
        for atom_idx in chain:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in carbon_adj[atom_idx]:
                if nbr in backbone_set:
                    continue
                # Validate branch from this neighbor (each attached branch is counted only once).
                valid, err = valid_branch(nbr, backbone_set)
                if not valid:
                    branch_errors.append(f"Invalid branch at backbone atom {atom_idx}: {err}")
                else:
                    branch_count += 1
        if branch_errors:
            continue  # chain candidate fails branch test.
        valid_candidates.append((chain, branch_count))
    
    if not valid_candidates:
        # If no candidate chain passed the branch test, try to supply one error reason from the longest overall chain.
        longest_chain = max(all_paths, key=len)
        if len(longest_chain) < 5:
            return False, "Main carbon chain is too short to be a fatty acid"
        else:
            return False, "No candidate main chain found that yields only methyl branches"
    
    # Choose the candidate with maximal main chain length.
    best_chain, best_branch_count = max(valid_candidates, key=lambda x: len(x[0]))
    main_chain_length = len(best_chain)
    msg = (f"CORRECT Methyl-branched fatty acid with {best_branch_count} methyl branch(es) "
           f"on a main chain of length {main_chain_length}")
    return True, msg

# Example usage: (if run as a script, one could test on a few SMILES)
if __name__ == '__main__':
    test_smiles_list = [
        "CC(C)CCCCCCCCCCCCCCCCC(O)=O",  # 18-methylnonadecanoic acid (expected true)
        "OC(=O)CCC(C)=C",              # 4-methyl-4-pentenoic acid (expected true)
        "CC(C)(C)C(O)=O",              # pivalic acid (expected false: chain too short)
        "OC(=O)CC(C(C)(C)C)=C",         # 3-tert-Butyl-3-butenoic acid (expected false if branch not allowed, but per our definition may be true)
    ]
    for smi in test_smiles_list:
        res, reason = is_methyl_branched_fatty_acid(smi)
        print(smi, res, reason)