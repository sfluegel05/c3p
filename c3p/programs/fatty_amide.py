"""
Classifies: CHEBI:29348 fatty amide
"""
#!/usr/bin/env python3
"""
Classifies: fatty amide
A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
Heuristic: presence of an amide functional group (C(=O)N) where the acyl side 
(i.e. the substituent attached to the carbonyl carbon that is not the amide N)
contains a contiguous carbon chain of a minimum length (here set as >=4 carbons).
"""

from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    
    In our heuristic we require the following:
      - The molecule must contain at least one amide moiety, identified by the substructure C(=O)N.
      - For at least one such amide group, the acyl substituent (attached to the carbonyl C apart from O and N)
        must contain a contiguous chain of carbon atoms of minimum length (>= 4 carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets our criteria for a fatty amide, False otherwise.
        str: Explanation of the criteria met or not met.
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an amide group: C(=O)N
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide (C(=O)N) functional group found"
    
    # Helper function: given an atom index to start from, traverse connected carbons only.
    def count_contiguous_carbons(start_idx, visited):
        count = 0
        stack = [start_idx]
        while stack:
            idx = stack.pop()
            if idx in visited:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            # Make sure the atom is carbon (atomic num 6)
            if atom.GetAtomicNum() != 6:
                continue
            count += 1
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in visited and nbr.GetAtomicNum() == 6:
                    stack.append(nbr_idx)
        return count

    # We now examine each amide match. We assume the order in the SMARTS "C(=O)N" is:
    # index 0: carbonyl carbon, index 1: oxygen, index 2: amide nitrogen.
    # For the carbonyl carbon, find neighbor(s) that are not the oxygen (index1) and not the amide N (index2).
    fatty_amide_found = False
    reason = ""
    for match in amide_matches:
        carbonyl_idx = match[0]
        oxygen_idx = match[1]
        amideN_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Get the neighboring atom indices of the carbonyl carbon.
        neighbors = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors()]
        # Exclude the oxygen and the amide nitrogen.
        acyl_candidates = [n for n in neighbors if n not in (oxygen_idx, amideN_idx)]
        
        if not acyl_candidates:
            # No acyl substituent found on the carbonyl carbon.
            reason = "Amide group present but no acyl substituent attached to the carbonyl carbon."
            continue
        
        # For each acyl candidate, count the number of contiguous carbon atoms in that substituent.
        # We want to check that the portion of the fatty chain (excluding the carbonyl carbon itself)
        # has at least 4 carbon atoms.
        for candidate in acyl_candidates:
            candidate_atom = mol.GetAtomWithIdx(candidate)
            # Only consider candidate if it is a carbon atom
            if candidate_atom.GetAtomicNum() != 6:
                continue
            # Perform a DFS over carbon atoms starting at the candidate.
            visited = set()
            chain_length = count_contiguous_carbons(candidate, visited)
            if chain_length >= 4:
                return True, f"Found amide group with an acyl chain of {chain_length} contiguous carbons."
            else:
                reason = f"Found amide group but acyl chain only has {chain_length} contiguous carbons (need at least 4)."
                
    return False, reason if reason else "No fatty acyl chain meeting criteria found"