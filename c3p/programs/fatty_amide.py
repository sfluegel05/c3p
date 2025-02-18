"""
Classifies: CHEBI:29348 fatty amide
"""
#!/usr/bin/env python3
"""
Classifies: fatty amide

A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
Heuristic improvements in this version:
  1. Identify an amide (C(=O)N) group.
  2. From the carbonyl carbon, exclude the oxygen and the amide nitrogen and choose the acyl candidate.
  3. Walk in a linear (non–branched) manner along eligible (aliphatic, nonaromatic, non‐ring) carbons.
     A qualifying acyl chain must contain at least 4 contiguous carbons.
  4. Verify that the acyl chain accounts for at least a minimum fraction (here 20%) of all heavy atoms.
  5. Check that the substituent attached to the amide N (besides the carbonyl partner) is “small” (fewer than 20 heavy atoms).
If no qualifying amide group is found, the program returns False with an explanation.
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines whether the given SMILES string corresponds to a fatty amide.
    The molecule must contain a C(=O)N group where:
      - The substituent on the carbonyl carbon (the acyl portion) is a contiguous linear chain
        of at least 4 aliphatic carbons (nonaromatic, non‐ring) and this chain represents a significant part of the molecule (at least 20% of heavy atoms).
      - Additionally, the substituent on the amide nitrogen (aside from the carbonyl carbon) should be relatively small (<20 heavy atoms).
      
    Args:
      smiles (str): The SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple of a boolean and an explanatory string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Build heavy atoms count for later: count atoms with atomic number > 1.
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)

    # Define a SMARTS to find an amide group: carbonyl carbon (index 0), carbonyl oxygen (index 1), and amide nitrogen (index 2)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide (C(=O)N) functional group found"
    
    # An iterative function to “walk” along a linear chain.
    def linear_chain_length(start_idx, parent_idx):
        """
        Walks (iteratively) from start_idx (attached to parent_idx) along a chain of carbon atoms.
        Only nonaromatic, nonring carbons are allowed. The walk stops if more than one eligible neighbor is encountered.
        Returns (length, chain_atom_indices) where length includes the starting atom.
        """
        chain = [start_idx]
        current_idx = start_idx
        prev_idx = parent_idx
        while True:
            current_atom = mol.GetAtomWithIdx(current_idx)
            # Get eligible neighbors: must be carbon, nonaromatic, not in a ring, and not the one we came from.
            eligible = []
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == prev_idx:
                    continue
                if nbr.GetAtomicNum() != 6:
                    continue
                if nbr.GetIsAromatic() or nbr.IsInRing():
                    continue
                eligible.append(nbr_idx)
            # For a linear (non‐branched) chain, there should be exactly one eligible neighbor.
            if len(eligible) == 1:
                chain.append(eligible[0])
                prev_idx, current_idx = current_idx, eligible[0]
            else:
                break
        return len(chain), chain

    # A DFS to count the size of the fragment (the amine substituent) attached to the amide nitrogen.
    def count_heavy_atoms(atom_idx, excluded, visited):
        visited = visited | {atom_idx}
        count = 1  # count this atom
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited or nbr_idx in excluded:
                continue
            if mol.GetAtomWithIdx(nbr_idx).GetAtomicNum() == 1:  # skip hydrogens
                continue
            count += count_heavy_atoms(nbr_idx, excluded, visited)
        return count

    # We set a threshold for the amine substituent size.
    AMINE_SIZE_THRESHOLD = 20
    # We also require that the acyl chain represents at least this fraction of all heavy atoms.
    ACYL_FRACTION_THRESHOLD = 0.20
    reasons = []
    
    for match in amide_matches:
        carbonyl_idx = match[0]
        oxy_idx = match[1]
        amideN_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Determine acyl candidate(s): the neighbors of the carbonyl carbon, except the oxygen and the amide nitrogen.
        neighbors = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors()]
        acyl_candidates = [n for n in neighbors if n not in (oxy_idx, amideN_idx)]
        if not acyl_candidates:
            reasons.append("Amide group present but no acyl substituent attached to the carbonyl carbon.")
            continue
        
        acyl_chain_ok = False
        acyl_reason = ""
        for candidate in acyl_candidates:
            candidate_atom = mol.GetAtomWithIdx(candidate)
            # Candidate must be carbon, nonaromatic, and not in a ring.
            if candidate_atom.GetAtomicNum() != 6:
                continue
            if candidate_atom.GetIsAromatic() or candidate_atom.IsInRing():
                continue
            # Walk along the supposed acyl chain.
            chain_length, chain_atoms = linear_chain_length(candidate, carbonyl_idx)
            if chain_length < 4:
                acyl_reason = f"Found amide group but acyl chain only has {chain_length} contiguous aliphatic carbons (need at least 4)."
                continue
            # Check that the chain represents a significant fraction of the molecule.
            if chain_length / total_heavy < ACYL_FRACTION_THRESHOLD:
                acyl_reason = f"Acyl chain ({chain_length} carbons) is too small a fraction of the molecule (fraction {chain_length/total_heavy:.2f}; needed ≥ {ACYL_FRACTION_THRESHOLD})."
                continue
            # If we pass both tests then we have a qualifying acyl chain.
            acyl_chain_ok = True
            acyl_reason = f"Found amide group with an acyl chain of {chain_length} contiguous aliphatic carbons."
            break
        
        if not acyl_chain_ok:
            reasons.append(acyl_reason)
            continue
        
        # Now check the N‐side (amine substituent) size.
        amideN_atom = mol.GetAtomWithIdx(amideN_idx)
        # Exclude the carbonyl carbon.
        n_neighbors = [nbr.GetIdx() for nbr in amideN_atom.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
        amine_size = 0
        for nbr in n_neighbors:
            amine_size += count_heavy_atoms(nbr, {carbonyl_idx}, set())
        if amine_size > AMINE_SIZE_THRESHOLD:
            reasons.append(f"Acyl chain qualifies, but the amine fragment attached to the amide nitrogen is too large ({amine_size} heavy atoms; threshold {AMINE_SIZE_THRESHOLD}).")
            continue
        
        # If we reach here for any amide match, we classify the molecule as a fatty amide.
        return True, acyl_reason

    if reasons:
        # Return the first encountered reason.
        return False, reasons[0]
    else:
        return False, "No fatty amide (with a qualifying fatty acyl chain) found."

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "NC(CCCCCCCCCCCCCCC)=O",         # hexadecanamide: expected True
        "O=C(NCCC1=CC=CC=C1)CCCC",         # N-(2-phenylethyl)pentanamide: expected True
        "CCCCCCCC(=O)NCCO",               # N-(octanoyl)ethanolamine: expected True
        "CCCCCCCCCCCC(N)=O",              # dodecanamide: expected True
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCCCC(=O)Nc1ccccc1",  # linolenic acid anilide: expected True
        "O=C(NCC[NH3+])C[C@](CC(=O)NC[C@@H](C([O-])=O)[NH3+])(C([O-])=O)O" # example false positive from peptides
    ]
    for smi in test_smiles:
        flag, reason = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {flag}\nReason: {reason}\n")