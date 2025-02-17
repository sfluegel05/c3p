"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: A medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any medium-chain fatty acid.

This program uses these heuristics:
  1. The molecule must contain a CoA moiety (we check via a characteristic SMARTS fragment).
  2. There must be a thioester group (C(=O)S) linking the fatty acyl chain to CoA.
  3. Starting from the carbonyl carbon (of the thioester) the longest simple carbon chain 
     (ignoring atoms in rings and excluding CoA fragment atoms) is computed.
  4. The chain (including the carbonyl carbon) must be between 6 and 12 carbons long.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a medium-chain fatty acyl-CoA.
    
    The definition requires that the molecule contains the CoA moiety and a thioester linkage
    that attaches a fatty acyl chain. The fatty acyl chain (starting at the carbonyl carbon and
    following only contiguous acyclic carbon (atomic number 6) atoms outside of the CoA fragment)
    must have a chain length between 6 and 12 carbons (inclusive).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as medium-chain fatty acyl-CoA, False otherwise.
        str: Explanation/reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Check for CoA moiety. We use a characteristic fragment SMARTS.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"
    
    # Get all atom indices that are part of any CoA fragment occurrence.
    coa_idxs = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        coa_idxs.update(match)
    
    # 2. Look for a thioester group: a carbonyl directly bonded to a sulfur.
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
        
    # Helper: Compute the longest linear (acyclic) chain of carbon atoms starting from a given atom.
    # Only traverse atoms with atomic number == 6 (carbon), that are not in any ring and not in the CoA fragment.
    def longest_chain_from(atom, visited):
        idx = atom.GetIdx()
        # If already visited, prevent cycles.
        if idx in visited:
            return 0
        # If atom is not carbon, or if it belongs to the CoA fragment, ignore.
        if atom.GetAtomicNum() != 6 or idx in coa_idxs:
            return 0
        # If atom is in a ring then do not continue (fatty acyl chains are acyclic).
        if atom.IsInRing():
            return 0
        
        visited.add(idx)
        max_length = 1  # count current atom.
        for nbr in atom.GetNeighbors():
            # Only consider neighbor carbons (and ensure we do not go back into CoA)
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in coa_idxs:
                length = 1 + longest_chain_from(nbr, visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    reasons = []
    valid_match_found = False
    
    # 3. For each found thioester, try to measure the fatty acyl chain length.
    for match in thioester_matches:
        # Per our SMARTS "C(=O)S": match[0] is the carbonyl carbon, match[1] is the sulfur.
        carbonyl_idx, sulfur_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # From the carbonyl carbon get all neighbors that are carbon and not the sulfur.
        acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                          if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum() == 6]
        if not acyl_neighbors:
            reasons.append("Thioester carbonyl does not have a neighboring alkyl carbon")
            continue
        
        # We assume the fatty acyl chain follows from one of these neighbors.
        # (In typical thioester, the carbonyl carbon has one carbon neighbor on the acyl side.)
        for acyl_start in acyl_neighbors:
            # Compute the longest acyclic chain starting from the carbonyl atom.
            # We count the carbonyl carbon itself plus the longest chain starting from the selected neighbor.
            chain_len = 1 + longest_chain_from(acyl_start, set())
            if 6 <= chain_len <= 12:
                reason = f"Found thioester linkage with fatty acyl chain length of {chain_len} carbons (medium-chain)."
                return True, reason
            else:
                reasons.append(f"Thioester found but fatty acyl chain length is {chain_len} carbons (not medium-chain).")
    
    # If none of the thioester linkages yielded a valid chain length:
    return False, " ; ".join(reasons)
    
# (End of program)