"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: A medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any medium-chain fatty acid.

Heuristics used:
 1. The molecule must contain a recognizable CoA moiety (identified by a characteristic SMARTS fragment).
 2. There must be a thioester group (C(=O)S) linking a fatty acyl chain to the CoA.
 3. Starting from the thioester carbonyl carbon (which is counted as carbon 1), the molecule is “walked”
    along a contiguous, acyclic chain. Only atoms that are carbon (atomic number 6) or sulfur (atomic number 16)
    (to capture thia–substituted analogs) are allowed, provided they are not part of any CoA fragment.
    Only bonds of type single or double are followed.
 4. The longest linear extension (if branches occur) is chosen.
 5. To be “medium‐chain”, the chain length (including the carbonyl) must be between 6 and 12.
  
Note that this is a heuristic and will mis‐classify some edge cases.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a medium-chain fatty acyl-CoA.
    
    It checks that:
       - a CoA moiety is present (via a characteristic substructure)
       - a thioester linkage (C(=O)S) exists
       - the acyl chain connected (starting from the thioester carbonyl carbon) extends 
         (following only linear, acyclic bonds through permitted atoms) for a length between 6 and 12 atoms.
         The allowed atoms in the chain are carbon (C) and, to catch thia replacements, sulfur (S).
         
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if classified as medium-chain fatty acyl-CoA, False otherwise.
       str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Check for CoA moiety.
    # Use a characteristic fragment that is commonly found in CoA.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"
    
    # Gather atom indices that are in any CoA fragment match.
    coa_idxs = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        coa_idxs.update(match)
    
    # 2. Look for a thioester linkage. We search for the C(=O)S pattern.
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
        
    # Helper: recursively determine the length of a contiguous linear chain.
    # We allow atoms that are carbon (6) or sulfur (16) as part of the acyl chain.
    # We only traverse bonds of type SINGLE or DOUBLE.
    def linear_chain_length(atom, prev_idx):
        length = 1  # count current atom
        next_extensions = []
        for bond in atom.GetBonds():
            # Ignore bonds not single or double
            if bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                continue
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetIdx() == prev_idx:
                continue
            # Only allow chain atoms that are C (6) or S (16), not in CoA, and not in a ring.
            if nbr.GetAtomicNum() in (6, 16) and (nbr.GetIdx() not in coa_idxs) and (not nbr.IsInRing()):
                next_extensions.append(linear_chain_length(nbr, atom.GetIdx()))
        if next_extensions:
            return length + max(next_extensions)
        else:
            return length

    reasons = []
    # 3. Process each thioester match.
    for match in thioester_matches:
        # Identify the carbonyl carbon in the thioester.
        # In the SMARTS "C(=O)S" the match tuple should include the carbonyl carbon;
        # we assume this is at index 0.
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Verify that one of its neighbors is sulfur (the thioester S).
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
        
        # The fatty acyl chain should be linked to the carbonyl carbon via a neighbor that is not:
        #   (a) the sulfur, and (b) part of the CoA fragment.
        acyl_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() 
                          if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum() in (6, 16)
                          and nbr.GetIdx() not in coa_idxs]
        if not acyl_neighbors:
            reasons.append("Thioester carbonyl does not have a neighboring alkyl (or thia-alkyl) atom outside CoA.")
            continue
        
        # Evaluate each potential acyl chain neighbor.
        # Count chain length starting from the carbonyl (which is counted as 1);
        # then add the length of the longest contiguous extension.
        for acyl_start in acyl_neighbors:
            chain_len = 1 + linear_chain_length(acyl_start, carbonyl_idx)
            if 6 <= chain_len <= 12:
                reason = f"Found thioester linkage with fatty acyl chain length of {chain_len} atoms (medium-chain)."
                return True, reason
            else:
                reasons.append(f"Thioester found but fatty acyl chain length is {chain_len} atoms (not between 6 and 12).")
    
    if reasons:
        return False, " ; ".join(reasons)
    return False, "No valid fatty acyl chain detected"

# (End of program)