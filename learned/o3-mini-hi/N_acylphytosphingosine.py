"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
#!/usr/bin/env python3
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
This program attempts to improve on the earlier heuristic by (1) explicitly identifying,
for a given amide bond (C(=O)N), the acyl branch (which must have a contiguous chain of at least 10 carbons)
and (2) the sphingoid (phytosphingosine-like) branch (which must show at least 6 contiguous carbons in its “tail”
plus evidence of two hydroxyl groups on the head region). The overall molecular weight is also required to be >500 Da.
Note: This heuristic approach may not cover all edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def longest_carbon_chain(mol, start_idx, exclude=set()):
    """
    Recursively computes the longest contiguous chain length (number of carbon atoms)
    starting from the atom with index start_idx. Only moves over carbon–carbon single bonds.
    """
    atom = mol.GetAtomWithIdx(start_idx)
    if atom.GetAtomicNum() != 6:
        return 0
    max_length = 1
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        nbr = bond.GetOtherAtom(atom)
        nbr_idx = nbr.GetIdx()
        if nbr_idx in exclude:
            continue
        if nbr.GetAtomicNum() == 6:
            new_exclude = set(exclude)
            new_exclude.add(start_idx)
            length = 1 + longest_carbon_chain(mol, nbr_idx, new_exclude)
            if length > max_length:
                max_length = length
    return max_length

def count_direct_hydroxyls(atom):
    """
    Counts oxygen atoms directly attached to the provided atom that are likely part of OH groups.
    (Checks that the bond is single and that the oxygen has at least one hydrogen.)
    """
    oh_count = 0
    for bond in atom.GetBonds():
        # Only consider single bonds
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        nbr = bond.GetOtherAtom(atom)
        if nbr.GetAtomicNum() == 8:
            # Count oxygen if it appears to be an –OH (by having at least one hydrogen)
            if nbr.GetTotalNumHs() >= 1:
                oh_count += 1
    return oh_count

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    Heuristic criteria applied for each amide bond (C(=O)N):
      1. From the carbonyl carbon, the acyl side must deliver a contiguous chain of >=10 carbons.
      2. From the amide nitrogen, the non-acyl branch is expected to be the sphingoid part.
         We require that its “tail” (longest contiguous carbon chain) be at least 6 carbons and that
         the head region (the alpha carbon plus one neighboring carbon) shows evidence of at least 2 -OH groups.
      3. The overall molecular weight is expected to be >500 Da.
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      bool: True if the molecule meets the criteria for N-acylphytosphingosine, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all general amide bonds using SMARTS "C(=O)N"
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found"
    
    # Loop over each amide match and try to classify properly.
    # In our SMARTS, we expect the match tuple to be (carbonyl C, carbonyl O, amide N)
    for match in amide_matches:
        carbonyl_C_idx, carbonyl_O_idx, amide_N_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_C_idx)
        amide_N_atom = mol.GetAtomWithIdx(amide_N_idx)
        
        # Identify acyl branch: among neighbors of the carbonyl carbon, choose a carbon
        # that is not the amide nitrogen or the carbonyl oxygen.
        acyl_branch_indices = []
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in (amide_N_idx, carbonyl_O_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_branch_indices.append(nbr_idx)
        if not acyl_branch_indices:
            continue  # no proper acyl branch found; try next amide
        
        # For simplicity, take the acyl branch with the longest chain.
        acyl_chain_length = 0
        for idx in acyl_branch_indices:
            chain_len = longest_carbon_chain(mol, idx, exclude={carbonyl_C_idx})
            if chain_len > acyl_chain_length:
                acyl_chain_length = chain_len
        if acyl_chain_length < 10:
            # Acyl branch too short for this amide; skip to next
            continue
        
        # Identify sphingoid branch: choose the neighbor of amide N that is carbon and not the carbonyl C.
        sph_branch = None
        for nbr in amide_N_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_C_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                sph_branch = nbr
                break
        if sph_branch is None:
            continue
        
        # Compute sphingoid tail length:
        sph_tail_length = longest_carbon_chain(mol, sph_branch.GetIdx(), exclude={amide_N_idx})
        if sph_tail_length < 6:
            continue
        
        # Now check for sufficient hydroxylation on the sphingoid head.
        # We check the alpha carbon (sph_branch) and one of its carbon neighbors (if present).
        oh_head_count = count_direct_hydroxyls(sph_branch)
        # Also check one other neighbor (if another carbon is attached) that is not the amide N.
        for nbr in sph_branch.GetNeighbors():
            if nbr.GetIdx() == amide_N_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                oh_head_count += count_direct_hydroxyls(nbr)
                break  # only consider one additional neighbor for the head region
        
        if oh_head_count < 2:
            continue
        
        # Check overall molecular weight (usually ceramides are >500 Da)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 500:
            return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a ceramide structure"
        
        # If we have reached here, we consider this amide bond as giving a valid N-acylphytosphingosine.
        reason = (f"Found amide bond connecting a long acyl chain (length {acyl_chain_length}) "
                  f"to a sphingoid-like unit with tail length {sph_tail_length} and "
                  f"{oh_head_count} hydroxyl(s) in the head region; molecular weight {mol_wt:.1f} Da")
        return True, reason
    
    # If none of the amide bonds satisfy all criteria, provide a rationale.
    return False, "No amide bond was found with a fatty acyl chain (>=10 C) and a sphingoid unit (>=6 C tail with >=2 OH) on the nitrogen"

# Example usage:
if __name__ == "__main__":
    # Test one provided example.
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
    valid, reason = is_N_acylphytosphingosine(test_smiles)
    print("Is N-acylphytosphingosine?", valid)
    print("Reason:", reason)