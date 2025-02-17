"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: A medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any medium-chain fatty acid.

Heuristics used:
 1. The molecule must contain a recognizable CoA moiety. We use a characteristic fragment 
    found in the CoA’s phosphoester part. (SMARTS: "COP(O)(=O)OP(O)(=O)")
 2. There must be a thioester linkage ("C(=O)S") joining the fatty acyl chain to CoA.
 3. From the thioester carbonyl carbon (counted as 1), a DFS is used to compute the length 
    of the contiguous acyclic chain (traversing only C or S atoms that are not part of CoA).
 4. For a medium‐chain fatty acyl, the computed chain length should be between 6 and 12.
    
Note: These heuristics may omit edge cases.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines whether the given molecule (as a SMILES string) is a 
    medium-chain fatty acyl-CoA according to heuristics.

    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a medium-chain fatty acyl-CoA, else False.
      str: Explanation (reason) for the result.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a CoA moiety by using a simplified and valid fragment SMARTS.
    coa_smarts = "COP(O)(=O)OP(O)(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in CoA SMARTS pattern"
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"
    
    # Get indices for atoms in any CoA fragment match.
    coa_idxs = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        coa_idxs.update(match)
        
    # 2. Look for a thioester linkage defined by the SMARTS "C(=O)S"
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Error in thioester SMARTS pattern"
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # 3. Helper DFS to find the longest simple acyclic chain.
    # Only traverse atoms that are allowed (atomic numbers {6 (C), 16 (S)}),
    # avoid atoms in the CoA moiety, and skip ring atoms.
    def dfs_longest_path(atom, coming_from, visited):
        longest = 0  # branch length count (excluding the start atom; caller adds 1)
        for bond in atom.GetBonds():
            # Only allow single/double bonds.
            if bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetIdx() == coming_from:
                continue
            if nbr.GetIdx() in visited:
                continue
            # Only traverse allowed atoms.
            if nbr.GetAtomicNum() not in (6, 16):
                continue
            # Do not traverse if the neighbor is part of the CoA moiety.
            if nbr.GetIdx() in coa_idxs:
                continue
            # Do not traverse ring atoms.
            if nbr.IsInRing():
                continue
            new_visited = visited.union({nbr.GetIdx()})
            branch_length = 1 + dfs_longest_path(nbr, atom.GetIdx(), new_visited)
            if branch_length > longest:
                longest = branch_length
        return longest
    
    reasons = []
    # 4. Process each thioester match.
    for match in thioester_matches:
        # Assume in the thioester SMARTS "C(=O)S", index 0 is the carbonyl carbon.
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Verify that the carbonyl has a neighboring sulfur.
        sulfur_found = False
        sulfur_idx = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetSymbol() == "S":
                sulfur_found = True
                sulfur_idx = nbr.GetIdx()
                break
        if not sulfur_found:
            reasons.append("Thioester match did not yield an adjacent sulfur.")
            continue
        
        # Identify neighbor atoms (outside the thioester S) that are allowed and not in CoA.
        acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                          if nbr.GetIdx() != sulfur_idx and 
                             nbr.GetAtomicNum() in (6, 16) and 
                             nbr.GetIdx() not in coa_idxs and 
                             not nbr.IsInRing()]
        if not acyl_neighbors:
            reasons.append("Thioester carbonyl has no valid acyl neighbor outside CoA.")
            continue
        
        for acyl_start in acyl_neighbors:
            # Start chain length counting with the carbonyl (1 atom)
            visited = {carbonyl_idx, acyl_start.GetIdx()}
            branch_length = 1 + dfs_longest_path(acyl_start, carbonyl_idx, visited)
            # Check if chain length (including carbonyl) is between 6 and 12.
            if 6 <= branch_length <= 12:
                reason = f"Found thioester linkage with fatty acyl chain length of {branch_length} atoms (medium-chain)."
                return True, reason
            else:
                reasons.append(f"Found thioester linkage but fatty acyl chain length is {branch_length} atoms (not between 6 and 12).")
    
    if reasons:
        return False, " ; ".join(reasons)
    return False, "No valid fatty acyl chain detected"

# End of program