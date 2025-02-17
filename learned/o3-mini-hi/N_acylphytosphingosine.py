"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
#!/usr/bin/env python3
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
Heuristic criteria for each candidate amide bond (C(=O)N):
  1. The acyl branch (from the carbonyl C) must have a contiguous carbon chain of at least 10 carbons.
  2. The sphingoid branch (from the amide N) must have a carbon tail of >=6 carbons and the head region must show at least 2 –OH groups.
  3. The molecular weight should be >500 Da.
Note: This approach is heuristic and may not cover every edge case.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol, start_idx, visited=None):
    """
    Recursively computes the longest contiguous carbon chain starting at atom index start_idx.
    Only traverses carbon–carbon single bonds.
    """
    if visited is None:
        visited = set()
    atom = mol.GetAtomWithIdx(start_idx)
    if atom.GetAtomicNum() != 6:  # only carbons are counted
        return 0
    visited.add(start_idx)
    max_length = 1  # count the starting carbon
    for bond in atom.GetBonds():
        # only consider single bonds between carbons
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        nbr = bond.GetOtherAtom(atom)
        nbr_idx = nbr.GetIdx()
        if nbr_idx in visited:
            continue
        if nbr.GetAtomicNum() == 6:
            chain_len = 1 + longest_carbon_chain(mol, nbr_idx, visited.copy())
            if chain_len > max_length:
                max_length = chain_len
    return max_length

def count_direct_hydroxyls(atom):
    """
    Counts the number of direct oxygen neighbors that appear to be part of an –OH group.
    (Checks that O is attached via a single bond and that it has at least one hydrogen.)
    """
    oh_count = 0
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        nbr = bond.GetOtherAtom(atom)
        if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
            oh_count += 1
    return oh_count

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    
    For each amide bond (recognized by SMARTS "C(=O)N") we:
      1. Identify the acyl branch starting from the carbonyl carbon (exclude the amide N and carbonyl O).
         The longest contiguous carbon chain must be >= 10 carbons.
      2. Identify the sphingoid branch from the amide nitrogen. For at least one carbon neighbor
         (other than the carbonyl carbon) the contiguous tail length must be >= 6 carbons.
         In addition, we count –OH groups on the candidate (the head region) and at one of its carbon neighbors.
         The sum should be >= 2.
      3. The overall molecular weight (exact MW) must be > 500 Da.
      
    Returns:
      (bool, str): (True, explanation) if criteria met, else (False, explanation) for failure.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all amide bonds via SMARTS "C(=O)N"
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found"

    # Loop over candidate amide bonds and try to classify.
    for match in amide_matches:
        # In our SMARTS, we expect: index0 = carbonyl C, index1 = carbonyl O, index2 = amide N
        carbonyl_C_idx, carbonyl_O_idx, amide_N_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_C_idx)
        amide_N_atom = mol.GetAtomWithIdx(amide_N_idx)

        # (1) Identify acyl branch on the carbonyl carbon:
        acyl_candidates = []
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Exclude the amide nitrogen and the carbonyl oxygen
            if nbr_idx in (amide_N_idx, carbonyl_O_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_candidates.append(nbr_idx)
        if not acyl_candidates:
            continue  # try next amide bond if no acyl branch

        # Find the maximum acyl chain length among candidates
        max_acyl_length = 0
        for idx in acyl_candidates:
            chain_len = longest_carbon_chain(mol, idx, visited={carbonyl_C_idx})
            if chain_len > max_acyl_length:
                max_acyl_length = chain_len
        if max_acyl_length < 10:
            # Acyl chain too short for this amide bond candidate.
            continue

        # (2) Identify sphingoid branch from the amide nitrogen:
        sph_candidates = []
        for nbr in amide_N_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_C_idx:
                continue  # skip the acyl side attachment
            if nbr.GetAtomicNum() == 6:
                sph_candidates.append(nbr)
        if not sph_candidates:
            continue  # no sphingoid branch candidate

        # For each sphingoid candidate, check tail length and hydroxylation on the "head region"
        valid_sph = False
        for candidate in sph_candidates:
            cand_idx = candidate.GetIdx()
            # Compute tail length starting from candidate; do not go back to the amide N.
            tail_length = longest_carbon_chain(mol, cand_idx, visited={amide_N_idx})
            if tail_length < 6:
                continue  # tail not long enough for this candidate

            # Count hydroxyl groups attached directly to the candidate (alpha position)
            oh_count = count_direct_hydroxyls(candidate)
            # Also consider one additional carbon neighbor (if exists) that is not coming from the amide N
            for nbr in candidate.GetNeighbors():
                if nbr.GetIdx() == amide_N_idx:
                    continue
                if nbr.GetAtomicNum() == 6:
                    oh_count += count_direct_hydroxyls(nbr)
                    break  # only consider one additional neighbor
            if oh_count < 2:
                continue  # insufficient hydroxylation in head region
            # This candidate qualifies as the sphingoid branch.
            valid_sph = True
            sph_tail_length = tail_length
            break

        if not valid_sph:
            continue  # try next amide bond

        # (3) Check overall molecular weight (>500 Da)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 500:
            return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a ceramide structure"

        reason = (f"Found amide bond with acyl chain length {max_acyl_length} and "
                  f"sphingoid branch tail length {sph_tail_length} with sufficient hydroxylation; "
                  f"molecular weight {mol_wt:.1f} Da")
        return True, reason

    # If none of the amide bonds satisfy all criteria:
    return False, ("No amide bond was found with a fatty acyl chain (>=10 C) "
                   "and a sphingoid unit (>=6 C tail with >=2 OH) on the nitrogen")

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
    valid, reason = is_N_acylphytosphingosine(test_smiles)
    print("Is N-acylphytosphingosine?", valid)
    print("Reason:", reason)