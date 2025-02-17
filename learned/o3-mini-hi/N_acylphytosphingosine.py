"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
Heuristic criteria for a candidate amide bond (C(=O)N):
  1. On the carbonyl carbon: the acyl branch (excluding the amide N and the carbonyl O) must have a linear carbon chain of >=10 carbons.
  2. On the amide nitrogen: one of the connected carbon atoms (sphingoid branch) must have a linear tail of >=6 carbons and at least 2 hydroxyl (–OH) groups 
     (counting those directly on the candidate plus on one of its neighboring carbons).
  3. Overall molecular weight must be >500 Da.

This revised version improves the recursive chain-length search by keeping track of visited atoms.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def linear_chain_length(mol, current_idx, visited=None):
    """
    Recursively computes the longest linear chain length starting at current_idx,
    following only single bonds between carbon atoms and avoiding atoms already visited.
    Returns the number of carbon atoms in the chain (including the starting atom).
    """
    atom = mol.GetAtomWithIdx(current_idx)
    if atom.GetAtomicNum() != 6:
        return 0  # Only count carbon atoms
    if visited is None:
        visited = set()
    # Copy visited for the current path and add the current index.
    visited = visited.union({current_idx})
    max_length = 1
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue  # only follow single bonds
        neighbor = bond.GetOtherAtom(atom)
        nid = neighbor.GetIdx()
        if nid in visited:
            continue
        if neighbor.GetAtomicNum() == 6:
            branch_length = 1 + linear_chain_length(mol, nid, visited)
            if branch_length > max_length:
                max_length = branch_length
    return max_length

def count_direct_hydroxyls(atom):
    """
    Counts the number of directly attached hydroxyl (–OH) groups to a given atom.
    Considers bonds to oxygen atoms that are bound to at least one hydrogen.
    """
    oh_count = 0
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        neighbor = bond.GetOtherAtom(atom)
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
            oh_count += 1
    return oh_count

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.

    For each amide bond (recognized by SMARTS "C(=O)N") the function performs:
      1. On the carbonyl carbon: among all connected atoms (excluding the amide N and the carbonyl O) the longest linear carbon chain must be >=10 carbons.
      2. On the amide nitrogen: among connected carbon atoms (excluding the carbonyl C),
         at least one candidate must have a tail (linear chain) of >=6 carbons and show sufficient –OH groups in the head region (>=2 OH groups counting those on the candidate and one adjacent carbon).
      3. Checks if the overall exact molecular weight is >500 Da.

    Returns:
        (bool, str): (True, explanation) if criteria are met, otherwise (False, explanation).
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find amide bonds using SMARTS "C(=O)N"
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found in the molecule"

    # Process each candidate amide bond match.
    for match in amide_matches:
        # The SMARTS expected indices: index0 = carbonyl C, index1 = carbonyl O, index2 = amide N.
        if len(match) != 3:
            continue  # Unexpected match, skip.
        carbonyl_C_idx, carbonyl_O_idx, amide_N_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_C_idx)
        amide_N_atom = mol.GetAtomWithIdx(amide_N_idx)

        # (1) Check the acyl side on the carbonyl carbon.
        acyl_candidates = []
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() in (amide_N_idx, carbonyl_O_idx):
                continue  # Skip amide N and carbonyl O.
            if nbr.GetAtomicNum() == 6:
                acyl_candidates.append(nbr.GetIdx())
        if not acyl_candidates:
            continue  # No candidate acyl branch for this amide, continue to next match.
        max_acyl_length = 0
        for idx in acyl_candidates:
            chain_len = linear_chain_length(mol, idx)
            if chain_len > max_acyl_length:
                max_acyl_length = chain_len
        if max_acyl_length < 10:
            # The acyl chain is too short; try next candidate amide bond.
            continue

        # (2) Check the sphingoid branch coming off the amide nitrogen.
        sphingoid_candidates = []
        for nbr in amide_N_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_C_idx:
                continue  # Skip the acyl branch already considered.
            if nbr.GetAtomicNum() == 6:
                sphingoid_candidates.append(nbr)
        if not sphingoid_candidates:
            continue  # No potential sphingoid branch.
        valid_sph = False
        sph_tail_length = 0
        for candidate in sphingoid_candidates:
            tail_length = linear_chain_length(mol, candidate.GetIdx())
            if tail_length < 6:
                continue  # Tail too short, try next candidate.
            # Check hydroxyl groups: count OH on the candidate and on one adjacent carbon (if any).
            oh_count = count_direct_hydroxyls(candidate)
            # Look at one carbon neighbor (excluding the amide N) for additional OH groups.
            for nbr in candidate.GetNeighbors():
                if nbr.GetIdx() == amide_N_idx:
                    continue
                if nbr.GetAtomicNum() == 6:
                    oh_count += count_direct_hydroxyls(nbr)
                    break  # Consider only one additional neighbor.
            if oh_count < 2:
                continue  # Insufficient OH groups.
            valid_sph = True
            sph_tail_length = tail_length
            break  # Found a valid sphingoid branch.
        if not valid_sph:
            continue  # No sphingoid candidate passed the criteria in this amide match.

        # (3) Check molecular weight.
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 500:
            return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a ceramide structure"
        
        reason = (f"Found amide bond with fatty acyl chain (length {max_acyl_length} C) and "
                  f"sphingoid branch (tail length {sph_tail_length} C with sufficient –OH groups); "
                  f"molecular weight {mol_wt:.1f} Da")
        return True, reason

    return False, ("No amide bond was found with a fatty acyl chain (>=10 C) "
                   "and a sphingoid unit (>=6 C tail with >=2 OH) on the nitrogen")

# Example usage:
if __name__ == "__main__":
    # Testing one example SMILES from the provided list:
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
    valid, reason = is_N_acylphytosphingosine(test_smiles)
    print("Is N-acylphytosphingosine?", valid)
    print("Reason:", reason)