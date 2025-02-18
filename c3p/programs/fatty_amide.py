"""
Classifies: CHEBI:29348 fatty amide
"""
#!/usr/bin/env python3
"""
Classifies: fatty amide

A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
In our approach we search for an amide group (C(=O)N) where:
  - One substituent on the carbonyl carbon (the acyl side) is a contiguous, linear,
    nonaromatic, non‐ring carbon chain of at least MIN_CHAIN_LENGTH carbons.
  - We compute a modified ratio:
         (acyl_chain_length) / (acyl_chain_length + aliphatic_count_on_N_side)
    where on the N side we only count contiguous aliphatic carbons.
  - The ratio must meet (or exceed) a threshold. This helps to select cases where
    the fatty (acyl) chain “dominates” the amide fragment.
If no qualifying amide group is found, the function returns (False, reason).
  
Note: By counting only nonaromatic, non‐ring carbon atoms on each side we discount
large peptide backbones and decorated N‐substituents that previously raised false positives.
"""

from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines whether the given SMILES string corresponds to a fatty amide.
    The molecule must contain an amide (C(=O)N) group where:
      - One substituent on the carbonyl carbon (the acyl side) is a contiguous,
        linear aliphatic chain (nonaromatic, non‐ring) of at least MIN_CHAIN_LENGTH carbons.
      - The ratio: acyl_chain_length / (acyl_chain_length + aliphatic_count_Nside)
        is at least RATIO_THRESHOLD.
      
    Args:
       smiles (str): The SMILES string of the molecule.
     
    Returns:
      (bool, str): A tuple (result, explanation).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define our criteria:
    MIN_CHAIN_LENGTH = 4  # minimum contiguous carbons (as in fatty acid R-group)
    RATIO_THRESHOLD = 0.50

    # First, find amide groups using a SMARTS pattern.
    # The pattern orders: carbonyl carbon, attached oxygen, and adjacent N.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide (C(=O)N) functional group found"
    
    # If multiple amide groups are present, that is a red flag for peptide bonds.
    if len(amide_matches) > 1:
        # We still check each one, but record a note if it appears peptide‐like.
        peptide_warning = True
    else:
        peptide_warning = False

    # Helper function: Follow a chain linearly (one continuation only)
    # starting from a candidate atom (attached to the carbonyl) if it is a carbon
    # that is nonaromatic and not in a ring.
    def linear_chain_length(start_idx, parent_idx):
        chain = [start_idx]
        current_idx = start_idx
        prev_idx = parent_idx
        while True:
            current_atom = mol.GetAtomWithIdx(current_idx)
            next_candidates = []
            # Examine neighbors except the atom we came from.
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == prev_idx:
                    continue
                if nbr.GetAtomicNum() != 6:
                    continue
                if nbr.GetIsAromatic() or nbr.IsInRing():
                    continue
                next_candidates.append(nbr_idx)
            # For a linear walk we continue only if exactly one candidate is found.
            if len(next_candidates) == 1:
                chain.append(next_candidates[0])
                prev_idx, current_idx = current_idx, next_candidates[0]
            else:
                break
        return len(chain), chain

    # Helper DFS: Count contiguous aliphatic carbon atoms (nonaromatic, not in ring)
    # starting from a given atom. This DFS traverses all connected such carbons.
    def dfs_aliphatic(atom_idx, excluded, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Only count if it is carbon, non-aromatic, not in ring.
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
            # Only follow if neighbor is carbon and satisfies aliphatic condition.
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                count += dfs_aliphatic(nbr_idx, excluded, visited)
        return count

    reasons = []  # collect reasons if no match qualifies

    # Now iterate over each found amide group and test it.
    for match in amide_matches:
        # In our SMARTS "C(=O)N":
        # match[0] is the carbonyl carbon, match[1] should be the oxygen,
        # and match[2] the nitrogen.
        carbonyl_idx = match[0]
        oxy_idx = match[1]
        amideN_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Look at neighbors on the carbonyl carbon that are not the oxygen or the amide nitrogen.
        neighbors = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors()]
        acyl_candidates = [n for n in neighbors if n not in (oxy_idx, amideN_idx)]
        if not acyl_candidates:
            reasons.append("Amide group present but no acyl substituent attached to the carbonyl carbon.")
            continue

        qualified_amide = False
        amide_reason = ""
        for candidate in acyl_candidates:
            candidate_atom = mol.GetAtomWithIdx(candidate)
            # The acyl chain must start with a carbon that is nonaromatic and not in a ring.
            if candidate_atom.GetAtomicNum() != 6:
                continue
            if candidate_atom.GetIsAromatic() or candidate_atom.IsInRing():
                continue

            chain_length, chain_atoms = linear_chain_length(candidate, carbonyl_idx)
            if chain_length < MIN_CHAIN_LENGTH:
                amide_reason = (f"Found amide group but acyl chain only has {chain_length} contiguous aliphatic carbons "
                                f"(need at least {MIN_CHAIN_LENGTH}).")
                continue

            # Now, for the amine substituent attached at the N,
            # we look at neighbors of the amide N ignoring the carbonyl carbon.
            amideN_atom = mol.GetAtomWithIdx(amideN_idx)
            n_neighbors = [nbr.GetIdx() for nbr in amideN_atom.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
            amine_aliphatic_count = 0
            visited_sum = set()
            for nbr_idx in n_neighbors:
                amine_aliphatic_count += dfs_aliphatic(nbr_idx, {carbonyl_idx}, visited_sum)
            # If no aliphatic atoms are found (for instance if the substituent is aromatic)
            # we set amine_aliphatic_count to 0.
            # Compute the ratio of acyl chain carbons to total aliphatic carbons involved.
            ratio = chain_length / (chain_length + amine_aliphatic_count) if (chain_length + amine_aliphatic_count) > 0 else 0

            # For reporting reasons, show ratio to two decimals.
            if ratio >= RATIO_THRESHOLD:
                amide_reason = (f"Found amide group with an acyl chain of {chain_length} contiguous aliphatic carbons "
                                f"and an acyl/total (aliphatic) substituent ratio of {ratio:.2f}.")
                qualified_amide = True
                break
            else:
                amide_reason = (f"Acyl chain of {chain_length} carbons found, but its ratio compared to the "
                                f"amine aliphatic fragment (size {amine_aliphatic_count}) is only {ratio:.2f} (needed ≥ {RATIO_THRESHOLD}).")
        if qualified_amide:
            # If the molecule shows additional amide groups (likely peptide bonds), add a note.
            note = " "
            if peptide_warning:
                note = "However, multiple amide groups were detected; check for peptide structure."
            return True, amide_reason + note
        else:
            reasons.append(amide_reason)
            # Continue checking other matches.
        
    if reasons:
        return False, reasons[0]
    else:
        return False, "No fatty amide (with a qualifying fatty acyl chain) found."

# Example usage: run only if executed as script.
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "NC(CCCCCCCCCCCCCCC)=O",         # hexadecanamide
        "O=C(NCCC1=CC=CC=C1)CCCC",         # N-(2-phenylethyl)pentanamide (expected now to be True)
        "CCCCCCCC(=O)NCCO",               # N-(octanoyl)ethanolamine
        "CCCCCCCCCCCC(N)=O",              # dodecanamide
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCCCC(=O)Nc1ccccc1",  # linolenic acid anilide
        "O=C(NCCCCN)CCCCCCCCCCCCCCC",     # Palmitoylputrescine
        "CCCCCC=CCC=CCCCCCCCC(=O)N",       # octadeca-9,12-dienamide
        "C(\CCCCCCCC(=O)N[C@H](C(O)=O)[C@H](CC)C)=C\CCCCCCCC",  # N-oleoyl-L-isoleucine
        # False positives (should be rejected):
        "O=C(N[C@@H](C(=O)O)CCC(=O)N)/C=C/C=C/C1=CC=C([N+]([O-])=NC2=CC=C(/C=C/C=C/C(=O)NC(C(=O)O)CCC(=O)N)C=C2)C=C1",  # Azoxymycin A
        "O=C(N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO)[C@@H](NC(=O)[C@@H](N)CCC(=O)N)CCCCN",  # Gln-Trp-Leu
        "CCCC(=O)Nc1ncnc2n(cnc12)[C@@H]1O[C@@H]2COP(O)(=O)O[C@H]2[C@H]1O",  # N(6)-butyryl-cAMP (should fail due to extremely short acyl chain)
        "S1C(=N[C@@H](C1)C=C)C[C@H]([C@H](C[C@@H](CCN(C(=O)[C@@H](CC)C)C)C)C)C",  # Kalkitoxin (expected False)
    ]
    for smi in test_smiles:
        flag, reason = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {flag}\nReason: {reason}\n")