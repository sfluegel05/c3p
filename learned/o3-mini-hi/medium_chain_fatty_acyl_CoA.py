"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: A medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any medium-chain fatty acid.

Heuristics used (revised):
 1. The molecule must contain a recognizable CoA moiety. Here we require a fragment that
    is typically present in CoA – we include a fragment of the pantetheine unit plus part of the adenosine.
 2. There must be a thioester linkage (C(=O)S) joining the fatty acyl chain to the CoA.
 3. From the thioester carbonyl carbon (which we count as position 1), we compute the longest
    simple acyclic path that goes only through “allowed” atoms (C or S only) and not going into 
    atoms that have been identified as part of the CoA moiety.
 4. To be “medium‐chain”, the chain length (including the carbonyl) must be between 6 and 12.
 
Note that these are heuristics and some edge cases will still be mis‐classified.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines whether the given molecule (as a SMILES string) is a medium-chain fatty acyl-CoA.
    
    The algorithm checks:
      - that a CoA moiety is present (by locating a characteristic CoA substructure)
      - that a thioester linkage (C(=O)S) exists
      - that one can find a contiguous acyclic chain (using only allowed atoms C and S)
        extending from the thioester carbonyl carbon with a total length (counting the carbonyl as 1)
        in the range 6 to 12.
        
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a medium-chain fatty acyl-CoA, else False.
      str: Explanation (reason) for the result.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for CoA moiety.
    # We use a SMARTS pattern aimed to capture a characteristic fragment from CoA.
    # We require a fragment containing parts of the pantetheine unit and an adenine ring.
    # (This is a compromise between coverage and avoiding spurious matches.)
    coa_smarts = "COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]"  # a fragment from the nucleotide-pantetheine unit
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"
    
    # Get indices for atoms in any CoA fragment match.
    coa_idxs = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        coa_idxs.update(match)
        
    # 2. Look for a thioester linkage defined by the SMARTS "C(=O)S"
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # 3. A helper DFS to find the longest simple acyclic path starting from a given atom.
    # We only traverse atoms that are allowed (atomic number in {6,16}) and not in the CoA fragment.
    # We also only follow bonds that are SINGLE or DOUBLE.
    def dfs_longest_path(atom, coming_from, visited):
        longest = 0  # count (not including the starting atom; caller adds 1)
        for bond in atom.GetBonds():
            # Only permitted bond types
            if bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetIdx() == coming_from:
                continue
            if nbr.GetIdx() in visited:
                continue
            # Only traverse allowed atoms: carbon (6) or sulfur (16)
            if nbr.GetAtomicNum() not in (6, 16):
                continue
            # Do not traverse if the atom is part of the CoA moiety
            if nbr.GetIdx() in coa_idxs:
                continue
            # For better confidence that we are following an acyclic path, do not visit ring atoms.
            if nbr.IsInRing():
                continue
            # Use DFS recursion.
            new_visited = visited.union({nbr.GetIdx()})
            branch_length = 1 + dfs_longest_path(nbr, atom.GetIdx(), new_visited)
            if branch_length > longest:
                longest = branch_length
        return longest
    
    reasons = []
    # 4. Process each thioester match and use the DFS to find the acyl chain length.
    for match in thioester_matches:
        # In the SMARTS "C(=O)S", assume index 0 is the carbonyl carbon.
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Verify that one neighbor is sulfur (the thioester S)
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
        
        # Identify neighbor atoms (other than S) that are allowed (part of acyl chain) and outside CoA.
        acyl_neighbors = [
            nbr for nbr in carbonyl_atom.GetNeighbors() 
            if nbr.GetIdx() != sulfur_idx 
               and nbr.GetAtomicNum() in (6, 16)
               and (nbr.GetIdx() not in coa_idxs)
               and (not nbr.IsInRing())
        ]
        if not acyl_neighbors:
            reasons.append("Thioester carbonyl does not have a neighboring alkyl (or thia-alkyl) atom outside CoA.")
            continue
        
        for acyl_start in acyl_neighbors:
            # The chain length is computed as 1 (for the carbonyl) plus the longest DFS path
            # starting from the acyl neighbor.
            visited = {carbonyl_idx, acyl_start.GetIdx()}
            branch_length = 1 + dfs_longest_path(acyl_start, carbonyl_idx, visited)
            # Check if the chain length (including carbonyl) falls between 6 and 12.
            if 6 <= branch_length <= 12:
                reason = f"Found thioester linkage with fatty acyl chain length of {branch_length} atoms (medium-chain)."
                return True, reason
            else:
                reasons.append(f"Thioester found but fatty acyl chain length is {branch_length} atoms (not between 6 and 12).")
    
    if reasons:
        return False, " ; ".join(reasons)
    return False, "No valid fatty acyl chain detected"

# (End of program)