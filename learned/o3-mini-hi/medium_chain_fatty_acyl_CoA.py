"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: A medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any medium-chain fatty acid.
Criteria (heuristic):
  1. The molecule must include a CoA moiety. We check for a characteristic CoA fragment.
  2. The molecule must have a thioester group (C(=O)S) linking a fatty acyl chain to the CoA.
  3. The fatty acyl chain is defined as the connected set of carbon atoms beginning from 
     the carbonyl carbon in the thioester (excluding the S side) and should be 6 to 12 C atoms long.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise.
        str: A reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Verify presence of a characteristic CoA moiety.
    # We use a fragment SMARTS that is common in CoA: "SCCNC(=O)CCNC(=O)"
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"

    # Look for thioester group: a carbonyl directly attached to a sulfur.
    # SMARTS: "C(=O)S"
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (fatty acyl) linkage found"

    # Helper function to count connected carbon atoms in the fatty acyl portion.
    def count_acyl_chain(carbon_atom, exclude_atoms):
        # DFS starting from 'carbon_atom': count only carbon atoms (atomic number 6)
        visited = set()
        stack = [carbon_atom]
        count = 0
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:
                count += 1
            else:
                continue  # Only traverse carbon atoms.
            # Traverse neighbors but skip any that are in exclude_atoms.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetIdx() not in exclude_atoms:
                    if nbr.GetAtomicNum() == 6:
                        stack.append(nbr)
        return count

    # For each thioester match, identify the acyl chain.
    # Each match returns a tuple of atom indices where:
    #   match[0] = carbonyl carbon, match[1] = sulfur.
    valid_match_found = False
    reasons = []
    for match in thioester_matches:
        carbonyl_idx, sulfur_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        # Get neighbors of the carbonyl carbon.
        neighbors = carbonyl_atom.GetNeighbors()
        # Exclude the sulfur atom and any oxygen(s) (the carbonyl O).
        acyl_neighbors = [nbr for nbr in neighbors if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum() == 6]
        if not acyl_neighbors:
            reasons.append("Thioester carbonyl does not have a neighboring alkyl carbon")
            continue  # Try next match if available.
        # Assume the first eligible neighbor represents the fatty acyl chain
        acyl_start = acyl_neighbors[0]
        # To be safe, exclude atoms that are part of the CoA fragment.
        # We already matched the CoA pattern. Get indices of CoA fragment.
        coa_matches = mol.GetSubstructMatches(coa_pattern)
        coa_idxs = set()
        for m in coa_matches:
            coa_idxs.update(m)
        # We also exclude the sulfur (already excluded).
        exclude_idxs = coa_idxs.union({sulfur_idx})
        # Count carbons in the connected acyl chain.
        # We also want to include the carbonyl carbon (it is part of the acyl group).
        acyl_chain_set = set()
        # Use a DFS that traverses only carbon atoms and avoids excluded atoms.
        def dfs(atom):
            if atom.GetIdx() in acyl_chain_set:
                return
            if atom.GetAtomicNum() != 6:
                return
            if atom.GetIdx() in exclude_idxs:
                return
            acyl_chain_set.add(atom.GetIdx())
            for nbr in atom.GetNeighbors():
                # Only traverse through carbon atoms (avoid S, O, etc.)
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in exclude_idxs:
                    dfs(nbr)
        # Start DFS from the carbonyl atom AND the chosen acyl neighbor.
        dfs(carbonyl_atom)
        dfs(acyl_start)
        chain_length = len(acyl_chain_set)
        # For medium-chain fatty acids, total carbons are usually between 6 and 12 (inclusive)
        if 6 <= chain_length <= 12:
            valid_match_found = True
            reason = f"Found thioester linkage with fatty acyl chain length of {chain_length} carbons (medium-chain)."
            return True, reason
        else:
            reasons.append(f"Thioester linkage found but fatty acyl chain length is {chain_length} carbons (not medium-chain).")
            
    # If none of the thioester groups were in the medium chain range:
    if not valid_match_found:
        return False, " ; ".join(reasons)
        
    # Fallback (should not be reached)
    return False, "Unable to classify as medium-chain fatty acyl-CoA"
    
# (End of program)