"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
#!/usr/bin/env python3
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
This improved function uses two strategies:
  (1) Look for a phytosphingosine-like substructure by a SMARTS that allows optional substituents 
      on the alpha carbon (to allow glycosylation).
  (2) If that fails, use a fallback check by inspecting the amide nitrogen's non-acyl branch 
      to see whether it displays at least one hydroxyl group on the first or second carbon 
      and that a contiguous aliphatic “tail” of at least 6 carbons is present.
Additional checks include that the fatty acyl chain (from the carbonyl) be long (>=10 carbons)
and that the overall molecule has a reasonable molecular weight (>500 Da).
Note: This is a heuristic approach and may not be perfect.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol, start_idx, exclude=set()):
    """
    Recursively computes the longest contiguous chain length (number of carbon atoms)
    starting from the atom with index start_idx. Only considers carbon–carbon bonds
    and avoids atoms in the 'exclude' set.
    """
    atom = mol.GetAtomWithIdx(start_idx)
    if atom.GetAtomicNum() != 6:
        return 0
    max_length = 1
    for nbr in atom.GetNeighbors():
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

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES.
    Heuristic criteria:
      1. The molecule must be valid.
      2. It must contain an amide bond (C(=O)N).
      3. The nitrogen of at least one amide bond is attached to a phytosphingosine‐like motif.
         We attempt first a SMARTS match of a core (allowing an arbitrary substituent on the
         alpha (first) carbon so that glycosylated variants are accepted). If that fails, we
         use a fallback: the branch from the amide nitrogen (non–acyl part) must show evidence
         of at least one hydroxyl (–OH) on the alpha carbon or its CH2 substituent and yield a
         contiguous aliphatic chain of at least 6 carbons (reflecting the sphingoid tail).
      4. The fatty acyl chain (attached to the carbonyl carbon) must be long (>=10 carbons).
      5. The overall molecular weight is expected to be >500 Da.
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      bool: True if the molecule meets the criteria for N-acylphytosphingosine, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find generic amide bonds using SMARTS "C(=O)N"
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found"
    
    # First, try to find the phytosphingosine-like backbone via SMARTS.
    # This pattern accepts a nitrogen attached to a carbon that has any substituent (to allow glycosylation)
    # and that carries an –OH, followed by two carbons each bearing an –OH.
    sphingo_smarts = Chem.MolFromSmarts("N[C]([*])(O)[C](O)[C](O)")
    sphingo_matches = mol.GetSubstructMatches(sphingo_smarts)
    
    selected_match = None  # will hold a tuple (amide_match, sphingo_match) if found

    # Loop over amide matches to see if the amide nitrogen is part of a sphingoid substructure.
    for amide in amide_matches:
        # In our amide SMARTS, we assume the tuple is (carbonyl C, carbonyl O, amide N).
        amide_n_idx = amide[2]
        # Check if any sphingo match has the first atom equal to this amide nitrogen.
        for sph in sphingo_matches:
            if sph[0] == amide_n_idx:
                selected_match = (amide, sph)
                break
        if selected_match is not None:
            break

    # If the SMARTS match did not work, do a fallback check.
    # For each amide, get the neighbor of the amide N (not the carbonyl carbon)
    if selected_match is None:
        fallback_found = False
        for amide in amide_matches:
            amide_n = mol.GetAtomWithIdx(amide[2])
            # Among the neighbors of the amide nitrogen, exclude the carbonyl carbon (amide[0])
            candidate_neighbors = [nbr for nbr in amide_n.GetNeighbors() if nbr.GetIdx() != amide[0]]
            for cand in candidate_neighbors:
                if cand.GetAtomicNum() != 6:
                    continue
                # Check for evidence of –OH in the immediate sphere
                oh_count = 0
                for subnbr in cand.GetNeighbors():
                    if subnbr.GetAtomicNum() == 8:
                        # Check that the oxygen is singly bonded (a typical –OH)
                        if subnbr.GetTotalNumHs() >= 1:
                            oh_count += 1
                # Also inspect a CH2 substituent if present (e.g. in a CH2OH group)
                for subnbr in cand.GetNeighbors():
                    if subnbr.GetAtomicNum() == 6 and subnbr.GetDegree() <= 4:
                        for subsub in subnbr.GetNeighbors():
                            if subsub.GetAtomicNum() == 8 and subsub.GetTotalNumHs() >= 1:
                                oh_count += 1
                # We require at least one (or two) hydroxyl groups in the head region.
                if oh_count >= 1:
                    # Accept this amide match; we take cand as the beginning of the sphingoid base.
                    selected_match = (amide, (amide_n.GetIdx(), cand.GetIdx()))
                    fallback_found = True
                    break
            if fallback_found:
                break

    if selected_match is None:
        return False, "No amide bond connecting a fatty acyl unit with a phytosphingosine-like (hydroxylated) unit found"
    
    # At this point we have an amide match and either a SMARTS sphingoid match or a fallback.
    amide_match, sph_match = selected_match
    # For the acyl part: from the carbonyl carbon (amide_match[0]) get its neighbor that is not the amide N or O.
    acyl_carbon_idx = amide_match[0]
    atom_carbonyl = mol.GetAtomWithIdx(acyl_carbon_idx)
    acyl_neighbors = []
    for nbr in atom_carbonyl.GetNeighbors():
        if nbr.GetIdx() not in (amide_match[1], amide_match[2]) and nbr.GetAtomicNum() == 6:
            acyl_neighbors.append(nbr.GetIdx())
    if not acyl_neighbors:
        return False, "No carbon substituent on the acyl carbon found"
    acyl_chain_length = 0
    for idx in acyl_neighbors:
        chain_len = longest_carbon_chain(mol, idx, exclude={acyl_carbon_idx})
        if chain_len > acyl_chain_length:
            acyl_chain_length = chain_len
    if acyl_chain_length < 10:
        return False, f"Acyl chain too short (found chain length {acyl_chain_length}, need at least 10 carbons)"
    
    # For the sphingoid (phytosphingosine) part:
    # If we got a SMARTS match, we expect the pattern to be (N, C1, C2, C3) and the sphingoid tail is presumed to extend 
    # from C3. Otherwise (fallback) we simply take the branch from the amide N that is not the acyl side.
    if len(sph_match) >= 4:
        sph_tail_start = sph_match[3]
    else:
        # fallback: use the candidate neighbor we picked (sph_match[1] from fallback)
        sph_tail_start = sph_match[1]
    # Find the contiguous chain length from the sph_tail_start. (This is our “tail” length.)
    sph_tail_length = longest_carbon_chain(mol, sph_tail_start, exclude={amide_match[2]})
    if sph_tail_length < 6:
        return False, f"Sphingoid tail too short (found chain length {sph_tail_length}, need at least 6 carbons)"
        
    # Check overall molecular weight (most ceramides are above ~500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a ceramide structure"
    
    # Passed all tests.
    return True, ("Found amide bond connecting a long acyl chain (length "
                  f"{acyl_chain_length}) to a phytosphingosine-like unit "
                  f"with tail length {sph_tail_length} and molecular weight {mol_wt:.1f} Da")
    
# Example usage:
if __name__ == "__main__":
    # Try one of the provided examples.
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
    valid, reason = is_N_acylphytosphingosine(test_smiles)
    print("Is N-acylphytosphingosine?", valid)
    print("Reason:", reason)