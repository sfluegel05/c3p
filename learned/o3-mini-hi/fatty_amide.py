"""
Classifies: CHEBI:29348 fatty amide
"""
#!/usr/bin/env python3
"""
Classifies: fatty amide

A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
We search for an amide (C(=O)N) group where:
  - One substituent attached to the carbonyl (acyl side) is a contiguous, linear,
    aliphatic (nonaromatic, non‐ring) chain of at least MIN_CHAIN_LENGTH carbons.
  - We then count the total aliphatic contiguous atoms (using a DFS) on the N‐side.
  - We compute the ratio: (acyl_chain_length) / (acyl_chain_length + Nsubstituent_aliphatic_count)
    and require that it is at least RATIO_THRESHOLD so that the fatty (acyl) chain “dominates”.
If no qualifying amide group is found, the function returns (False, reason).
Note: if multiple amide groups are detected, a warning is appended that it might be
due to peptide-like structure.
"""

from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines whether the given SMILES string corresponds to a fatty amide.

    The molecule must contain an amide (C(=O)N) group where:
      - The acyl substituent attached to the carbonyl carbon is a contiguous, unbranched,
        nonaromatic, non‐ring aliphatic chain of at least MIN_CHAIN_LENGTH atoms.
      - The ratio: acyl_chain_length / (acyl_chain_length + N_side_aliphatic_count)
        is at least RATIO_THRESHOLD.
     
    Args:
       smiles (str): SMILES string of the molecule.
     
    Returns:
      (bool, str): Tuple with a boolean result and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Define our criteria:
    MIN_CHAIN_LENGTH = 4    # must have at least 4 contiguous aliphatic carbons on the acyl side
    RATIO_THRESHOLD = 0.50   # required ratio of (acyl chain) / (acyl chain + N-side chain)

    # --- Find amide groups using a SMARTS pattern.
    # The SMARTS "C(=O)N" picks the carbonyl carbon, its oxygen, and the N atom.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide (C(=O)N) functional group found"
    
    # If multiple amide groups are present, that might point to peptides.
    peptide_warning = (len(amide_matches) > 1)

    # --- Helper: walk along a chain from a starting atom (child) in a linear (unbranched) fashion.
    # We continue only if exactly one eligible (nonaromatic, non‐ring, carbon) neighbor is found.
    def linear_chain_length(start_idx, parent_idx):
        chain_length = 1  # count the current atom
        current_idx = start_idx
        prev_idx = parent_idx
        while True:
            current_atom = mol.GetAtomWithIdx(current_idx)
            # Find eligible neighbor(s) that are carbons, nonaromatic, not in a ring,
            # and that are not the atom we came from.
            next_candidates = []
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == prev_idx:
                    continue
                if nbr.GetAtomicNum() != 6:
                    continue
                if nbr.GetIsAromatic() or nbr.IsInRing():
                    continue
                next_candidates.append(nbr_idx)
            # For a linear chain, accept only if exactly one candidate continues the chain.
            if len(next_candidates) == 1:
                chain_length += 1
                prev_idx, current_idx = current_idx, next_candidates[0]
            else:
                break
        return chain_length

    # --- Helper DFS: count all contiguous aliphatic carbons (nonaromatic, non‐ring)
    # starting from a given atom. This DFS will follow all branches.
    def dfs_aliphatic(atom_idx, excluded, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.IsInRing():
            return 0
        if atom_idx in visited:
            return 0
        visited.add(atom_idx)
        count = 1  # count this atom
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in excluded or nbr_idx in visited:
                continue
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                count += dfs_aliphatic(nbr_idx, excluded, visited)
        return count

    reasons = []  # if no match qualifies, collect reasons

    # --- Process each detected amide group.
    for match in amide_matches:
        # According to our SMARTS "C(=O)N":
        # match[0] = carbonyl carbon; match[1] = oxygen; match[2] = amide nitrogen.
        carbonyl_idx = match[0]
        oxy_idx = match[1]
        amideN_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Identify substituents on the carbonyl carbon that are not the oxygen or amide N.
        neighbors = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors()]
        acyl_candidates = [n for n in neighbors if n not in (oxy_idx, amideN_idx)]
        if not acyl_candidates:
            reasons.append("Amide group present but no acyl substituent attached to the carbonyl carbon.")
            continue

        qualified_amide = False
        amide_reason = ""
        # --- For each acyl candidate, check if it meets the linear, nonbranched criteria.
        for candidate in acyl_candidates:
            candidate_atom = mol.GetAtomWithIdx(candidate)
            # Ensure the candidate is a carbon and is not aromatic/inring.
            if candidate_atom.GetAtomicNum() != 6 or candidate_atom.GetIsAromatic() or candidate_atom.IsInRing():
                continue
            # Walk along the chain from this candidate.
            chain_length = linear_chain_length(candidate, carbonyl_idx)
            if chain_length < MIN_CHAIN_LENGTH:
                amide_reason = (f"Found amide group but acyl chain only has {chain_length} contiguous aliphatic carbons "
                                f"(need at least {MIN_CHAIN_LENGTH}).")
                continue

            # --- Count contiguous aliphatic carbons on the N‐side (amine substituent).
            amideN_atom = mol.GetAtomWithIdx(amideN_idx)
            n_neighbors = [nbr.GetIdx() for nbr in amideN_atom.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
            n_aliphatic_count = 0
            visited = set()
            for nbr_idx in n_neighbors:
                n_aliphatic_count += dfs_aliphatic(nbr_idx, {carbonyl_idx}, visited)
            # Compute the ratio: (acyl chain carbons) / (acyl chain + N‐side aliphatic carbons)
            total = chain_length + n_aliphatic_count
            ratio = chain_length / total if total > 0 else 0
            # Report the ratio rounded to two decimals.
            if ratio >= RATIO_THRESHOLD:
                amide_reason = (f"Found amide group with an acyl chain of {chain_length} contiguous aliphatic carbons "
                                f"and an acyl/total (aliphatic) substituent ratio of {ratio:.2f}.")
                qualified_amide = True
                break
            else:
                amide_reason = (f"Acyl chain of {chain_length} carbons found, but its ratio compared to the "
                                f"amine aliphatic fragment (size {n_aliphatic_count}) is only {ratio:.2f} (needed ≥ {RATIO_THRESHOLD}).")
        if qualified_amide:
            note = ""
            if peptide_warning:
                note = " However, multiple amide groups were detected; check for peptide structure."
            return True, amide_reason + note
        else:
            reasons.append(amide_reason)
    
    if reasons:
        return False, reasons[0]
    else:
        return False, "No fatty amide (with a qualifying fatty acyl chain) found."

# Example usage: run only when executed as script.
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "NC(CCCCCCCCCCCCCCC)=O",         # hexadecanamide
        "O=C(NCCC1=CC=CC=C1)CCCC",         # N-(2-phenylethyl)pentanamide
        "CCCCCCCC(=O)NCCO",               # N-(octanoyl)ethanolamine
        "CCCCCCCCCCCC(N)=O",              # dodecanamide
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCCCC(=O)Nc1ccccc1",  # linolenic acid anilide
        "O=C(NCCCCN)CCCCCCCCCCCCCCC",     # Palmitoylputrescine
        "CCCCCC=CCC=CCCCCCCCC(=O)N",       # octadeca-9,12-dienamide
        "C(\\CCCCCCCC(=O)N[C@H](C(O)=O)[C@H](CC)C)=C\\CCCCCCCC",  # N-oleoyl-L-isoleucine
        # False positives (should be rejected):
        "O=C(N[C@@H](C(=O)O)CCC(=O)N)/C=C/C=C/C1=CC=C([N+]([O-])=NC2=CC=C(/C=C/C=C/C(=O)NC(C(=O)O)CCC(=O)N)C=C2)C=C1",  # Azoxymycin A
        "O=C(N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO)[C@@H](NC(=O)[C@@H](N)CCC(=O)N)CCCCN",  # Gln-Trp-Leu
        # False negatives (should be rejected due to short acyl chain or low ratio):
        "CCCC(=O)Nc1ncnc2n(cnc12)[C@@H]1O[C@@H]2COP(O)(=O)O[C@H]2[C@H]1O",  # N(6)-butyryl-cAMP: acyl chain too short (3 Cs)
        "S1C(=N[C@@H](C1)C=C)C[C@H]([C@H](C[C@@H](CCN(C(=O)[C@@H](CC)C)C)C)C)C",  # Kalkitoxin: acyl chain only 1 C contiguous
    ]
    for smi in test_smiles:
        result, explanation = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")