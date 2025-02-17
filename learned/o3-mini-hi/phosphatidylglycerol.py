"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: Phosphatidylglycerol 
Definition: A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.
Criteria (improved):
  1. The molecule must contain exactly one phosphorus atom.
  2. The phosphorus atom should be bound (via single bonds) to three oxygen atoms.
     However, we ignore any oxygen that is a dead end (i.e. bonded only to P).
  3. Among the remaining oxygen branches:
       • Exactly one branch should contain no acyl ester groups and show glycerol-like connectivity.
         (This is the phosphoglycerol headgroup branch.)
       • Exactly one branch should contain exactly two acyl ester groups (i.e. “OC(=O)” fragments)
         that are fully contained in that branch. (This is the diacylglycerol branch.)
If any of these conditions are not met, an explanatory reason is returned.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    
    The algorithm:
      1. Parse the molecule and ensure exactly one phosphorus atom is present.
      2. For the phosphorus atom, scan its single-bond oxygen neighbors.
         Ignore oxygen atoms which are “dead ends” (degree==1) as they are not part of a branch.
      3. For each candidate branch (found by doing a breadth-first search from that oxygen
         and not going back to the phosphorus), count the number of acyl ester groups,
         defined by the pattern "OC(=O)" and check for a glycerol-like connectivity.
      4. Expect one branch with exactly 2 acyl ester groups (the diacylglycerol branch)
         and one branch with 0 acyl esters where at least one carbon appears to be bonded
         to two oxygens (the glycerol headgroup branch).
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a phosphatidylglycerol.
      str: Explanation of the classification decision.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check exactly one phosphorus atom (atomic number 15)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom, found {len(p_atoms)}"
    p_atom = p_atoms[0]
    
    # 2. Get oxygen neighbors attached by a single bond.
    candidate_oxygen_branches = []
    for bond in p_atom.GetBonds():
        # We only consider single bonds (ignore the P=O double bond)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        neighbor = bond.GetOtherAtom(p_atom)
        # Only consider oxygen atoms
        if neighbor.GetAtomicNum() != 8:
            continue
        # Ignore oxygen if it is a dead end (only connected to phosphorus)
        if neighbor.GetDegree() <= 1:
            continue
        candidate_oxygen_branches.append(neighbor)
    
    # For a typical phosphatidylglycerol we expect exactly two branches:
    # one diacylglycerol branch and one headgroup branch.
    if len(candidate_oxygen_branches) != 2:
        return False, f"Expected 2 nontrivial oxygen branches from phosphorus, found {len(candidate_oxygen_branches)}"

    # Helper: Given a starting atom, get the set of atom indices reachable (branch)
    # without traversing a set of forbidden atoms (here, the phosphorus).
    def bfs_branch(start_atom, forbidden_idxs):
        branch = set()
        stack = [start_atom]
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in branch:
                continue
            branch.add(atom.GetIdx())
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in forbidden_idxs:
                    continue
                # Always traverse single bonds; we do not restrict due to bond order.
                stack.append(nb)
        return branch

    forbidden = {p_atom.GetIdx()}
    
    # Define ester SMARTS: a fragment "OC(=O)". We will check that the oxygen in the match lies
    # within the branch.
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    
    # Helper: Count acyl ester groups within a set of atom indices (branch)
    def count_acyl_esters(branch_idxs):
        count = 0
        matches = mol.GetSubstructMatches(ester_pattern)
        for match in matches:
            # The match tuple is (O, C, O); we count the ester only if the first O is in our branch.
            if match[0] in branch_idxs:
                count += 1
        return count

    # Helper: Check if the branch contains at least one carbon with at least two oxygen neighbors
    # (a very loose indicator of a glycerol substructure).
    def looks_like_glycerol(branch_idxs):
        for idx in branch_idxs:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # carbon
                oxy_count = 0
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() in branch_idxs:
                        oxy_count += 1
                if oxy_count >= 2:
                    return True
        return False

    headgroup_branch_found = False
    acyl_branch_found = False
    headgroup_branch_reason = ""
    acyl_branch_reason = ""
    
    # Process each candidate branch:
    for o_atom in candidate_oxygen_branches:
        branch_idxs = bfs_branch(o_atom, forbidden)
        ester_count = count_acyl_esters(branch_idxs)
        # If branch contains acyl esters, we expect exactly 2 for the diacylglycerol branch.
        if ester_count == 2:
            if acyl_branch_found:
                return False, "More than one branch with 2 acyl ester groups found"
            acyl_branch_found = True
            acyl_branch_reason = f"Found diacyl branch with {ester_count} acyl ester groups."
        # If branch lacks acyl ester groups, try to check for glycerol connectivity.
        elif ester_count == 0:
            if looks_like_glycerol(branch_idxs):
                if headgroup_branch_found:
                    return False, "More than one branch appears to be a glycerol headgroup"
                headgroup_branch_found = True
                headgroup_branch_reason = "Found glycerol headgroup branch (no acyl ester groups and glycerol-like connectivity)."
            else:
                # If branch has 0 esters but does not look like glycerol, record information.
                headgroup_branch_reason += "A branch with 0 acyl esters did not show glycerol-like connectivity. "
        else:
            # If a branch has a number of esters not equal to 0 or 2, we do not expect that.
            return False, f"Unexpected acyl ester count ({ester_count}) in one branch."
    
    if not headgroup_branch_found:
        return False, f"Phosphoglycerol headgroup branch not conclusively identified. {headgroup_branch_reason}"
    if not acyl_branch_found:
        return False, f"Diacylglycerol branch (with exactly 2 acyl ester groups) not found. {acyl_branch_reason}"
    
    return True, ("Molecule contains one phosphorus atom with a unique phosphoglycerol headgroup branch "
                  "and a diacylglycerol branch featuring exactly 2 acyl ester groups.")

# Example usage (testing with PG(8:0/8:0)):
if __name__ == "__main__":
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCC)COC(=O)CCCCCCC)(OC[C@@H](O)CO)(O)=O"
    result, reason = is_phosphatidylglycerol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)