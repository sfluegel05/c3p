"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
Heuristic criteria for a candidate amide bond (C(=O)N):
  1. The acyl branch on the carbonyl carbon (excluding the amide N and the carbonyl O) must have a linear carbon chain of >=10 carbons.
  2. The sphingoid branch on the amide nitrogen must have a carbon tail (linear chain) of >=6 carbons.
     Additionally, the “head region” should exhibit at least 2 hydroxyl (–OH) groups (counting those on the candidate carbon and one additional neighboring carbon).
  3. The overall molecular weight should be >500 Da.

This approach is heuristic and does not cover every edge case.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Helper function that computes the longest linear chain (only along single C-C bonds) ignoring the branch we came from.
def linear_chain_length(mol, start_idx, parent_idx):
    """
    Recursively computes the longest linear chain length beginning at start_idx,
    disallowing a backward step to parent_idx.
    Returns the number of carbon atoms in the chain (including the starting atom).
    """
    atom = mol.GetAtomWithIdx(start_idx)
    # Only count carbon atoms.
    if atom.GetAtomicNum() != 6:
        return 0
    max_len = 1
    for bond in atom.GetBonds():
        # Only consider single bonds.
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        neighbor = bond.GetOtherAtom(atom)
        if neighbor.GetIdx() == parent_idx:
            continue
        if neighbor.GetAtomicNum() == 6:
            branch_length = 1 + linear_chain_length(mol, neighbor.GetIdx(), start_idx)
            if branch_length > max_len:
                max_len = branch_length
    return max_len

def count_direct_hydroxyls(atom):
    """
    Counts the number of directly attached hydroxyl (–OH) groups to the given atom.
    Checks single bonds to oxygen atoms that have at least one hydrogen.
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
      1. On the carbonyl carbon: among all connected atoms (excluding the amide N and the carbonyl O),
         the longest linear carbon chain must be >=10 carbons (fatty acyl chain).
      2. On the amide nitrogen: among carbon neighbors (other than the carbonyl C),
         at least one candidate must have a tail length (linear chain) of >=6 carbons and show
         sufficient –OH groups in the head region (>=2 when counting direct –OH’s on the candidate
         plus on one additional carbon neighbor).
      3. Checks if the overall exact molecular weight is >500 Da.
    
    Returns:
      (bool, str): (True, explanation) if criteria are met; otherwise (False, explanation).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all amide bonds using a SMARTS query "C(=O)N"
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found in the molecule"

    # Loop over all candidate amide bonds.
    for match in amide_matches:
        # In our SMARTS match we expect indices: 
        #  index0 = carbonyl C, index1 = carbonyl O, index2 = amide N.
        carbonyl_C_idx, carbonyl_O_idx, amide_N_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_C_idx)
        amide_N_atom = mol.GetAtomWithIdx(amide_N_idx)

        # (1) Determine acyl branch on the carbonyl carbon.
        acyl_candidates = []
        for nbr in carbonyl_atom.GetNeighbors():
            # Exclude the connected amide nitrogen and the carbonyl oxygen.
            if nbr.GetIdx() in (amide_N_idx, carbonyl_O_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_candidates.append(nbr.GetIdx())
        if not acyl_candidates:
            continue  # Try the next amide bond if no acyl branch exists.

        max_acyl_length = 0
        for idx in acyl_candidates:
            chain_len = linear_chain_length(mol, idx, carbonyl_C_idx)
            if chain_len > max_acyl_length:
                max_acyl_length = chain_len
        if max_acyl_length < 10:
            # This acyl chain is too short.
            continue

        # (2) Determine sphingoid branch from the amide nitrogen.
        sphingoid_candidates = []
        for nbr in amide_N_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_C_idx:
                continue  # Skip the acyl side.
            if nbr.GetAtomicNum() == 6:
                sphingoid_candidates.append(nbr)
        if not sphingoid_candidates:
            continue  # No candidate sphingoid branch.

        valid_sph = False
        sph_tail_length = 0
        for candidate in sphingoid_candidates:
            # Compute the tail (linear) chain length starting from candidate (excluding amide N).
            tail_length = linear_chain_length(mol, candidate.GetIdx(), amide_N_idx)
            if tail_length < 6:
                continue
            # Count –OH groups on this candidate and, if available, on one of its branch carbon neighbors.
            oh_count = count_direct_hydroxyls(candidate)
            for nbr in candidate.GetNeighbors():
                if nbr.GetIdx() == amide_N_idx:
                    continue
                if nbr.GetAtomicNum() == 6:
                    oh_count += count_direct_hydroxyls(nbr)
                    break  # Only consider one adjacent carbon.
            if oh_count < 2:
                continue
            valid_sph = True
            sph_tail_length = tail_length
            break

        if not valid_sph:
            continue  # Sphingoid branch conditions were not met.

        # (3) Check the overall molecular weight requirement.
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
    # Test with one example from the provided list:
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
    valid, reason = is_N_acylphytosphingosine(test_smiles)
    print("Is N-acylphytosphingosine?", valid)
    print("Reason:", reason)