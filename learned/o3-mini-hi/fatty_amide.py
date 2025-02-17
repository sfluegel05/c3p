"""
Classifies: CHEBI:29348 fatty amide
"""
#!/usr/bin/env python3
"""
Classifies: fatty amide

A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
Heuristic improvements in this version:
  1. Identify an amide (C(=O)N) group.
  2. For the carbonyl carbon, select the substituent (acyl chain) that is not the carbonyl oxygen or the amide nitrogen.
  3. Verify that starting from that acyl carbon there exists a contiguous linear (non–branched) chain of at least 4 carbon atoms.
     Here “contiguous” means we only traverse carbon atoms that are not in rings and not aromatic.
  4. Also check that the amine substituent (the group attached to the amide nitrogen apart from the carbonyl carbon)
     is not too large (i.e. has fewer than a fixed threshold of heavy atoms), because fatty amides are typically derived
     from a fatty acid (large acyl chain) and a relatively small amine.
If no amide group meeting these criteria is found, the program returns False with an explanation.
"""

from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines whether the given SMILES string corresponds to a fatty amide.
    
    The molecule must contain an amide group, C(=O)N, where the substituent on the carbonyl side
    (i.e. the acyl chain) is a contiguous linear (noncyclic, nonaromatic) chain of at least 4 carbon atoms.
    Additionally, we require that the substituent on the amide nitrogen is not very large (here < 20 heavy atoms)
    so that the molecule is more likely to be derived from a fatty acid.
    
    Args:
      smiles (str): The SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple with True and a reason if classified as fatty amide, otherwise False and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an amide SMARTS; note that this returns a match where:
    # index0 = carbonyl C, index1 = carbonyl O, index2 = amide N.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide (C(=O)N) functional group found"

    # Recursive DFS to find the longest simple (non-repeating) path starting from a given atom.
    # We only allow traversal over carbon atoms that are non-aromatic and not in any ring.
    def longest_linear_path(atom_idx, visited):
        # Mark the current atom as visited.
        visited = visited | {atom_idx}
        max_len = 1  # count current atom
        atom = mol.GetAtomWithIdx(atom_idx)
        # Explore only neighbors that are carbon, non‐aromatic, and not in a ring.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIsAromatic() or nbr.IsInRing():
                continue
            # Recursively compute path length from neighbor.
            path_len = 1 + longest_linear_path(nbr_idx, visited)
            if path_len > max_len:
                max_len = path_len
        return max_len

    # A helper to count the heavy atoms (non‐hydrogen) in the fragment connected to a given starting atom,
    # while excluding a given atom (or set of atoms). This DFS does not discriminate by element.
    def count_heavy_atoms(atom_idx, excluded, visited):
        visited = visited | {atom_idx}
        count = 1  # count this atom
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited or nbr_idx in excluded:
                continue
            if nbr.GetAtomicNum() == 1:  # skip hydrogens
                continue
            count += count_heavy_atoms(nbr_idx, excluded, visited)
        return count

    # We set a threshold for the amine substituent size (on the N-side)
    AMINE_SIZE_THRESHOLD = 20
    reasons = []
    # Check each amide match.
    for match in amide_matches:
        carbonyl_idx = match[0]
        oxy_idx = match[1]
        amideN_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # From the carbonyl carbon, get “acyl” candidate(s) excluding the oxygen and the amide nitrogen.
        neighbors = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors()]
        acyl_candidates = [n for n in neighbors if n not in (oxy_idx, amideN_idx)]
        if not acyl_candidates:
            reasons.append("Amide group present but no acyl substituent attached to the carbonyl carbon.")
            continue
        acyl_chain_ok = False
        acyl_reason = ""
        for candidate in acyl_candidates:
            candidate_atom = mol.GetAtomWithIdx(candidate)
            if candidate_atom.GetAtomicNum() != 6:
                continue
            # Only consider candidate if it is aliphatic: not aromatic and not in ring.
            if candidate_atom.GetIsAromatic() or candidate_atom.IsInRing():
                continue
            chain_length = longest_linear_path(candidate, set())
            if chain_length >= 4:
                acyl_chain_ok = True
                acyl_reason = f"Found amide group with an acyl chain of {chain_length} contiguous aliphatic carbons."
                break
            else:
                acyl_reason = f"Found amide group but acyl chain only has {chain_length} contiguous aliphatic carbons (need at least 4)."
        if not acyl_chain_ok:
            reasons.append(acyl_reason)
            continue

        # Now, check the size of the substituent on the amide nitrogen (the N-side) to ensure it is not too bulky.
        amideN_atom = mol.GetAtomWithIdx(amideN_idx)
        # Get the atoms attached to the nitrogen aside from the carbonyl carbon.
        n_neighbors = [nbr.GetIdx() for nbr in amideN_atom.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
        if not n_neighbors:
            amine_size = 0
        else:
            # Count heavy atoms in the connected fragment beginning from the first neighbor.
            # (If there are multiple disconnected fragments we sum the counts separately.)
            amine_visited = set()
            amine_size = 0
            for nbr in n_neighbors:
                if nbr in amine_visited:
                    continue
                amine_size += count_heavy_atoms(nbr, {carbonyl_idx}, set())
                # We do not update amine_visited precisely, but the threshold is heuristic.
        if amine_size > AMINE_SIZE_THRESHOLD:
            reasons.append(f"Acyl chain qualifies, but the amine fragment ({amine_size} heavy atoms) is too large (threshold {AMINE_SIZE_THRESHOLD}).")
            continue

        # If we reach here, we assume we've found a fatty amide.
        return True, acyl_reason

    # If no match met all criteria, return False with an explanation.
    if reasons:
        return False, reasons[0]
    else:
        return False, "No fatty amide (with a qualifying fatty acyl chain) found."

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "NC(CCCCCCCCCCCCCCC)=O",         # hexadecanamide
        "O=C(NCCC1=CC=CC=C1)CCCC",         # N-(2-phenylethyl)pentanamide
        "CCCCCCCC(=O)NCCO",               # N-(octanoyl)ethanolamine
        "CCCCCCCCCCCC(N)=O",              # dodecanamide
        "CC(=O)NCCC1=CC(=C(C=C1)O)O"       # Example of N-acetyldopamine (should fail; acyl chain too short)
    ]
    for smi in test_smiles:
        flag, reason = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {flag}\nReason: {reason}\n")