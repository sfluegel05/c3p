"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: A medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any medium-chain fatty acid.

This program uses these heuristics:
  1. The molecule must contain a CoA moiety (detected via a characteristic SMARTS fragment).
  2. There must be a thioester group (C(=O)S) linking the fatty acyl chain to CoA.
  3. Starting from the thioester carbonyl carbon (ignoring the sulfur from CoA), 
     we "walk" along a linear, acyclic chain of carbon atoms 
     (ignoring atoms that belong to CoA) and count how many carbons are in the chain.
  4. The chain length (counting the carbonyl carbon as 1) must be between 6 and 12, inclusive.
  
If the chain contains a branch, we choose the longest linear extension.
Note that this is a heuristic and may still mis‐classify some edge cases.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a medium-chain fatty acyl-CoA.
    
    The definition requires that the molecule contains the CoA moiety and a thioester linkage(s)
    that attaches a fatty acyl chain. We then walk from the carbonyl carbon (of the thioester)
    along a contiguous chain (only carbons, not in rings and not in the CoA fragment) and require
    that the chain length (including the carbonyl carbon) is between 6 and 12.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as medium-chain fatty acyl-CoA, False otherwise.
      str: Explanation/reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Check for CoA moiety.
    # We use a characteristic fragment SMARTS for a part of CoA.
    # (This pattern is not unique but should capture a sizeable part of CoA.)
    coa_smarts = "SCCNC(=O)CCNC(=O)" 
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"
    
    # Get all atom indices that occur in any match to the CoA fragment.
    coa_idxs = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        coa_idxs.update(match)
        
    # 2. Look for a thioester linkage: carbonyl (C(=O)) directly bonded to sulfur.
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
        
    # Helper function: starting from a given carbon atom, walk along a linear (acyclic) chain.
    # The walk only follows carbon (atomic number 6) atoms that are not in rings and
    # are not part of the CoA moiety. The walk stops if a branch (multiple unvisited valid neighbors)
    # is reached; in that case we take the maximum extension found.
    def linear_chain_length(atom, prev_idx):
        # We count the current atom (which should be carbon)
        length = 1
        next_candidates = []
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == prev_idx:
                continue
            # Only follow carbon atoms that are not in a ring and not in CoA fragment.
            if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() not in coa_idxs) and (not nbr.IsInRing()):
                next_candidates.append(nbr)
        # If no next candidate, then chain ends.
        if not next_candidates:
            return length
        else:
            # In most fatty acyl chains there is only one way forward.
            branch_lengths = []
            for nbr in next_candidates:
                branch_lengths.append(linear_chain_length(nbr, atom.GetIdx()))
            return length + max(branch_lengths)
    
    reasons = []
    
    # 3. Process each thioester match
    for match in thioester_matches:
        # The SMARTS "C(=O)S" returns a match tuple where:
        # match[0] is the carbonyl carbon, match[1] is the oxygen of the carbonyl (implicitly)
        # match[?] is the sulfur; however, the SMARTS order may vary.
        # To be safe, we search the neighbors of the carbonyl carbon to find the sulfur.
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
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
        
        # The fatty acyl chain should be attached to the carbonyl carbon *not via the sulfur*.
        # We exclude the sulfur neighbor and any neighbor that is part of CoA.
        acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                          if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum()==6 and nbr.GetIdx() not in coa_idxs]
        if not acyl_neighbors:
            reasons.append("Thioester carbonyl does not have a neighboring alkyl carbon (outside CoA).")
            continue
        
        # In typical fatty-acyl chains, there is one neighbor leading away from the thioester.
        for acyl_start in acyl_neighbors:
            # Count chain length starting from the carbonyl (which we count as 1)
            chain_len = 1 + linear_chain_length(acyl_start, carbonyl_idx)
            if 6 <= chain_len <= 12:
                reason = f"Found thioester linkage with fatty acyl chain length of {chain_len} carbons (medium-chain)."
                return True, reason
            else:
                reasons.append(f"Thioester found but fatty acyl chain length is {chain_len} carbons (not in 6-12 range).")
    
    if reasons:
        return False, " ; ".join(reasons)
    return False, "No valid fatty acyl chain detected"

# (End of program)