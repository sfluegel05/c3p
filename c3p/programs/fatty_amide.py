"""
Classifies: CHEBI:29348 fatty amide
"""
#!/usr/bin/env python3
"""
Classifies: fatty amide

A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
In this version we look for an amide group (C(=O)N) with two substituents:
  - One side (the acyl chain) must be a contiguous, linear, non‐aromatic, non‐ring chain of at least 4 carbons.
  - The acyl chain must “dominate” the amide fragment, meaning that the ratio:
        acyl_chain_size / (acyl_chain_size + amine_substituent_size)
    is at least 0.50.
This aims both to ensure that a genuine fatty acid residue is present and to avoid false positives
from peptide bonds or heavily decorated amine substituents.
If no qualifying amide group is found, the function returns (False, reason).
"""

from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines whether the given SMILES string corresponds to a fatty amide.
    The molecule must contain an amide (C(=O)N) group where:
      - One substituent on the carbonyl carbon (the acyl side) is a contiguous,
        linear aliphatic chain (nonaromatic, non‐ring) of at least 4 carbon atoms.
      - The ratio of acyl chain heavy atoms to total (acyl chain + N‐substituent heavy atoms)
        is at least 0.50.
    
    Args:
       smiles (str): The SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple (result, explanation).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count overall heavy atoms (Z>1) if needed later.
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    
    # Define a simple SMARTS to find an amide group: carbonyl carbon, its oxygen, and an adjacent nitrogen.
    # The match orders are: carbonyl carbon, oxygen, nitrogen.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide (C(=O)N) functional group found"
    
    # Helper function: Walk linearly from a starting carbon (attached to the carbonyl carbon)
    # along a chain that is nonaromatic, non‐ring and composed solely of carbon atoms.
    def linear_chain_length(start_idx, parent_idx):
        """
        Walks the chain from start_idx (which is attached to parent_idx) as far as possible in a linear fashion.
        The allowed atoms are carbons that are not aromatic and not in rings.
        It stops if branching (more than one eligible neighbor) is encountered.
        Returns (chain_length, chain_indices)
        """
        chain = [start_idx]
        current_idx = start_idx
        prev_idx = parent_idx
        while True:
            current_atom = mol.GetAtomWithIdx(current_idx)
            eligible = []
            # Check all neighbors except the one we came from.
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == prev_idx:
                    continue
                # The candidate must be a carbon; allow both sp2 and sp3 provided not aromatic and not in a ring.
                if nbr.GetAtomicNum() != 6:
                    continue
                if nbr.GetIsAromatic() or nbr.IsInRing():
                    continue
                eligible.append(nbr_idx)
            # For a linear chain, ensure there is at most one eligible continuation.
            if len(eligible) == 1:
                chain.append(eligible[0])
                prev_idx, current_idx = current_idx, eligible[0]
            else:
                break
        return len(chain), chain

    # Helper function: recursively count heavy atoms (atomic number > 1) starting from atom_idx.
    def dfs_count(atom_idx, excluded, visited):
        visited.add(atom_idx)
        count = 1  # count self
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited or nbr_idx in excluded:
                continue
            if nbr.GetAtomicNum() == 1:  # skip hydrogens
                continue
            count += dfs_count(nbr_idx, excluded, visited)
        return count

    # We set our threshold for the ratio:
    RATIO_THRESHOLD = 0.50
    # And minimum acyl chain contiguous carbon count:
    MIN_CHAIN_LENGTH = 4
    
    reasons = []  # collect reasons for failing different matches

    # For each amide match, try to find a qualifying acyl substituent.
    for match in amide_matches:
        carbonyl_idx = match[0]
        oxy_idx = match[1]
        amideN_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Determine acyl candidate(s): neighbors of the carbonyl carbon that are not oxygen or amide nitrogen.
        neighbors = [nbr.GetIdx() for nbr in carbonyl_atom.GetNeighbors()]
        acyl_candidates = [n for n in neighbors if n not in (oxy_idx, amideN_idx)]
        if not acyl_candidates:
            reasons.append("Amide group present but no acyl substituent attached to the carbonyl carbon.")
            continue
        
        acyl_chain_ok = False
        acyl_reason = ""
        for candidate in acyl_candidates:
            candidate_atom = mol.GetAtomWithIdx(candidate)
            # Must be a carbon, not aromatic, and not in a ring.
            if candidate_atom.GetAtomicNum() != 6:
                continue
            if candidate_atom.GetIsAromatic() or candidate_atom.IsInRing():
                continue
            # Walk the chain starting from this candidate.
            chain_length, chain_atoms = linear_chain_length(candidate, carbonyl_idx)
            if chain_length < MIN_CHAIN_LENGTH:
                acyl_reason = (f"Found amide group but acyl chain only has {chain_length} contiguous aliphatic carbons "
                               f"(need at least {MIN_CHAIN_LENGTH}).")
                continue
            
            # Now compute the amine substituent’s heavy atom count.
            # Start from the amide nitrogen and follow all neighbors except the carbonyl.
            amideN_atom = mol.GetAtomWithIdx(amideN_idx)
            n_neighbors = [nbr.GetIdx() for nbr in amideN_atom.GetNeighbors() if nbr.GetIdx() != carbonyl_idx]
            amine_size = 0
            visited_overall = set()
            for nbr in n_neighbors:
                amine_size += dfs_count(nbr, {carbonyl_idx}, visited_overall.copy())
            # (If no substituent exists, consider it as a very small fragment.)
            if amine_size == 0:
                amine_size = 1
            # Compute the ratio of the acyl chain heavy atoms to the combined substituents (acyl + amine).
            ratio = chain_length / (chain_length + amine_size)
            # Accept if the ratio is at least our threshold.
            if ratio >= RATIO_THRESHOLD:
                acyl_chain_ok = True
                acyl_reason = (f"Found amide group with an acyl chain of {chain_length} contiguous aliphatic carbons "
                               f"and an acyl/total substituent ratio of {ratio:.2f}.")
                break
            else:
                acyl_reason = (f"Acyl chain of {chain_length} carbons found, but its ratio compared to the "
                               f"amine fragment (size {amine_size}) is only {ratio:.2f} (needed ≥ {RATIO_THRESHOLD}).")
        if not acyl_chain_ok:
            reasons.append(acyl_reason)
            continue
        
        # If we reach here for any amide match, we classify the molecule as a fatty amide.
        return True, acyl_reason
    
    if reasons:
        return False, reasons[0]
    else:
        return False, "No fatty amide (with a qualifying fatty acyl chain) found."

# Example usage: (When run as a script, try a few test cases.)
if __name__ == '__main__':
    test_smiles = [
        "NC(CCCCCCCCCCCCCCC)=O",         # hexadecanamide: expected True
        "O=C(NCCC1=CC=CC=C1)CCCC",         # N-(2-phenylethyl)pentanamide: expected True
        "CCCCCCCC(=O)NCCO",               # N-(octanoyl)ethanolamine: expected True
        "CCCCCCCCCCCC(N)=O",              # dodecanamide: expected True
        "CC\\C=C/C\\C=C/C\\C=C/CCCCCCCC(=O)Nc1ccccc1",  # linolenic acid anilide: expected True
        # Some examples that were previously false positives/negatives:
        "O=C(N[C@H](C(=O)N[C@@H](C(C)C)C(O)=O)CO)[C@@H](N)CCCCN", # Lys-Ser-Val (false positive previously)
        "CCCC(=O)Nc1ncnc2n(cnc12)[C@@H]1O[C@@H]2COP(O)(=O)O[C@H]2[C@H]1O", # N(6)-butyryl-cAMP (expected False due to short acyl chain)
        "N1C=C(CCNC(CCCCCCC/C=C\\CCCCCCCC)=O)C2=C1C=CC(O)=C2", # N-oleoylserotonin (expected True)
        "S1C(=N[C@@H](C1)C=C)C[C@H]([C@H](C[C@@H](CCN(C(=O)[C@@H](CC)C)C)C)C)C", # Kalkitoxin (expected False)
    ]
    for smi in test_smiles:
        flag, reason = is_fatty_amide(smi)
        print(f"SMILES: {smi}\nResult: {flag}\nReason: {reason}\n")