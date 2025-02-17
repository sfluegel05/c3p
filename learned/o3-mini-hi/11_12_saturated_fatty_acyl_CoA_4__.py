"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
CHEBI:84947

This classifier verifies that the molecule is a fatty acyl-CoA(4-)
and that its fatty acyl portion – numbered linearly starting with the thioester
carbonyl as C1 – has a saturated (single) bond between C11 and C12.

Improvements over the previous version:
  • Instead of simply taking the longest DFS path, we enumerate all simple carbon-only
    paths emerging from the thioester carbonyl (after excluding the sulfur).
  • For each candidate chain we compute a branch penalty (extra carbon “branches” along the path)
    so that the best candidate is that which best corresponds to a linear acyl chain.
  • If multiple candidate chains disagree on the bond type between C11 and C12 then we err on
    the side of caution.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-).

    The molecule must contain a fatty acyl-CoA(4-) component. First, we search
    for a proxy substructure for CoA and a thioester linkage. Then, starting from the
    thioester carbonyl we perform a DFS (only following carbon neighbors) to enumerate
    all possible linear acyl chain paths. To decide which chain is the physiologically
    relevant fatty acyl chain, we assign a penalty to each candidate chain based on any extra
    branching relative to a fully linear chain. (A linear chain interior carbon should have exactly
    two carbon neighbors on the chain; additional carbon neighbors are regarded as branches.)
    We then choose the candidate with the lowest penalty (and if tied, the longest chain).
    Finally, numbering the complete fatty acyl chain (thioester carbonyl = C1, then the chosen carbon chain
    as C2, C3, etc), we require that the bond between C11 and C12 (i.e. chain atoms at indices 10 and 11)
    is a single (saturated) bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple where the first element is True if the molecule meets the criteria,
                     and the second element is a reason message.
    """
    # Parse SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a CoA fragment using a proxy SMARTS.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found – not a fatty acyl-CoA(4-) molecule"
    
    # Look for a thioester linkage: carbonyl ([CX3](=O)) directly attached to a sulfur ([SX2]).
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester group not found – molecule may not be fatty acyl-CoA(4-)"
    
    # Assume the first thioester match gives the fatty acyl linkage.
    # In our SMARTS the first atom is the carbonyl, the second is the sulfur.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the neighbor of the carbonyl atom that is (a) not the sulfur and (b) a carbon.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() == sulfur_idx:
            continue
        if nbr.GetAtomicNum() == 6:
            acyl_start = nbr
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found at the thioester group"
    
    # We now want to enumerate all simple acyl chain paths starting from acyl_start. We follow only carbon atoms.
    # To avoid revisiting atoms, we use a DFS that tracks visited atom indices.
    # To limit computation, we also set a maximum path length.
    MAX_DEPTH = 50
    def dfs(current_atom, visited, depth):
        # Returns a list of paths; each path is a list of atom objects starting with current_atom.
        paths = [[current_atom]]
        if depth >= MAX_DEPTH:
            return paths
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited:
                continue
            new_visited = visited | {nbr.GetIdx()}
            sub_paths = dfs(nbr, new_visited, depth+1)
            for sp in sub_paths:
                paths.append([current_atom] + sp)
        return paths

    all_paths = dfs(acyl_start, {acyl_start.GetIdx()}, 1)
    if not all_paths:
        return False, "No acyl chain path found"
    
    # Now, for each DFS path, we define the complete fatty acyl chain as:
    #  [carbonyl_atom] + path.
    # We then compute a "branch penalty" for the carbon-only chain defined as follows:
    # For each internal atom (not first or last), count how many carbon neighbors (in the entire molecule)
    # exist other than the two that are the immediate neighbors in the chain.
    # A perfectly linear chain should yield a penalty of 0.
    def chain_penalty(full_chain):
        penalty = 0
        # full_chain[0] is the carbonyl; full_chain[1:] are the chain carbons.
        for i in range(1, len(full_chain) - 1):
            atom = full_chain[i]
            # Identify neighbors that are carbons.
            cnbrs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            # The chain connects this atom to a previous and next atom (if present).
            expected = 0
            if i-1 >= 0:
                expected += 1
            if i+1 < len(full_chain):
                expected += 1
            extra = len(cnbrs) - expected
            if extra > 0:
                penalty += extra
        return penalty

    candidate_chains = []
    for path in all_paths:
        # The full chain always begins with the thioester carbonyl.
        full_chain = [carbonyl_atom] + path
        if len(full_chain) >= 2:
            pen = chain_penalty(full_chain)
            candidate_chains.append((full_chain, pen))
            
    if not candidate_chains:
        return False, "No candidate fatty acyl chain found"
    
    # Choose candidate(s) with lowest penalty. If more than one, choose the longest among these.
    min_pen = min(pen for (_, pen) in candidate_chains)
    best_candidates = [chain for (chain, pen) in candidate_chains if pen == min_pen]
    best_chain = max(best_candidates, key=lambda c: len(c))
    
    # For a fatty acyl chain to have an 11-12 bond, it must have at least 12 atoms (including carbonyl as C1).
    if len(best_chain) < 12:
        return False, f"Fatty acyl chain is too short ({len(best_chain)} atoms including the carbonyl) to have an 11-12 bond"
    
    # Because of branch alternatives, one may get more than one candidate chain with the same minimal penalty.
    # We now check all candidate chains with penalty == min_pen and maximal length.
    max_len = max(len(c) for (c, p) in candidate_chains if p==min_pen)
    final_candidates = [c for (c, p) in candidate_chains if p==min_pen and len(c)==max_len]
    
    # For each candidate, check the bond type between C11 (index 10) and C12 (index 11).
    # (Note: indexing – carbonyl is C1, then acyl carbons follow.)
    bond_types = set()
    for chain in final_candidates:
        atom_11 = chain[10]  # C11
        atom_12 = chain[11]  # C12
        bond = mol.GetBondBetweenAtoms(atom_11.GetIdx(), atom_12.GetIdx())
        if bond is None:
            bond_types.add("None")
        else:
            bond_types.add(bond.GetBondType())
    
    # If all candidate chains have a single bond between C11 and C12, we classify it as 11,12-saturated.
    if all(bt == Chem.BondType.SINGLE for bt in bond_types):
        return True, "The fatty acyl chain has a saturated (single) bond between carbon 11 and 12"
    else:
        return False, f"At least one candidate acyl chain shows a non-single (unsaturated) bond between C11 and C12: {bond_types}"

# Example usage: testing with one true positive SMILES.
if __name__ == "__main__":
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, ":", reason)