"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
#!/usr/bin/env python3
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
This function attempts to verify that the molecule contains an amide bond (C(=O)N) where the N is part
of a phytosphingosine‐like fragment (N[C@@H](CO)[C@H](O)[C@H](O)…) and that both the acyl chain and the sphingoid tail 
are reasonably long.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol, start_idx, exclude=set()):
    """
    Recursively computes the longest contiguous chain length (number of carbons) 
    starting from the atom with index start_idx. Only considers carbon atoms (atomic number 6)
    and does a simple DFS avoiding atoms in the 'exclude' set.
    """
    atom = mol.GetAtomWithIdx(start_idx)
    # Only count if carbon
    if atom.GetAtomicNum() != 6:
        return 0
    max_length = 1  # count the starting carbon
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx in exclude:
            continue
        if nbr.GetAtomicNum() == 6:
            # Mark current as visited along this branch
            new_exclude = set(exclude)
            new_exclude.add(start_idx)
            length = 1 + longest_carbon_chain(mol, nbr_idx, new_exclude)
            if length > max_length:
                max_length = length
    return max_length

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    Criteria (heuristic):
      1. The molecule must be valid.
      2. It must contain at least one amide bond (C(=O)N).
      3. The nitrogen of the amide bond must be embedded in a sphingoid-like motif 
         (heuristic match for a phytosphingosine backbone, e.g., N[C@@H](CO)[C@H](O)[C@H](O)).
      4. The acyl (fatty acid) part attached to the carbonyl carbon should be a long chain.
      5. The sphingoid (long-chain base) “tail” attached to the sphingoid core should be long.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule meets the criteria for N-acylphytosphingosine, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # SMARTS for a generic amide: C(=O)N
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide bond (C(=O)N) found"
        
    # SMARTS for a phytosphingosine-like core.
    # We ignore stereochemistry in this substructure (the symbols @@ and @) as they are not mandatory in the input.
    sphingo_smarts = Chem.MolFromSmarts("N[C](CO)[C](O)[C](O)")
    sphingo_matches = mol.GetSubstructMatches(sphingo_smarts)
    if not sphingo_matches:
        return False, "No phytosphingosine backbone (N[C](CO)[C](O)[C](O)) found"
        
    # Now look for an amide match in which the N atom overlaps with the sphingoid backbone.
    matching_found = False
    selected_amide = None
    selected_sphingo = None
    for amide in amide_matches:
        # amide tuple is (carbonyl C, oxygen, N)
        amide_n = amide[2]
        for sphingo in sphingo_matches:
            # sphingo tuple is (N, C1, C2, C3) where the first atom is the N of the sphingoid base.
            if sphingo[0] == amide_n:
                matching_found = True
                selected_amide = amide
                selected_sphingo = sphingo
                break
        if matching_found:
            break
    if not matching_found:
        return False, "No amide bond connecting a fatty acyl moiety with a phytosphingosine backbone found"
    
    # Now, check that the fatty acyl (acyl chain) part is long.
    # The acyl part is the portion attached to the carbonyl carbon (amide match index 0).
    acyl_carbon = selected_amide[0]
    # From the carbonyl carbon, the neighbors include the oxygen and the amide N. 
    # Choose the neighbor that is not the oxygen or the N.
    acyl_neighbors = []
    atom_carbonyl = mol.GetAtomWithIdx(acyl_carbon)
    for nbr in atom_carbonyl.GetNeighbors():
        if nbr.GetIdx() not in (selected_amide[1], selected_amide[2]):
            if nbr.GetAtomicNum() == 6:  # only consider carbons
                acyl_neighbors.append(nbr.GetIdx())
    if not acyl_neighbors:
        return False, "No carbon substituent on the acyl carbon found"
    
    # For each candidate neighbor, compute the length of the contiguous carbon chain.
    acyl_chain_length = 0
    for nbr_idx in acyl_neighbors:
        chain_len = longest_carbon_chain(mol, nbr_idx, exclude={acyl_carbon})
        if chain_len > acyl_chain_length:
            acyl_chain_length = chain_len
    # Set a minimum threshold (e.g. at least 10 carbons) for the acyl chain.
    if acyl_chain_length < 10:
        return False, f"Acyl chain too short (found chain length {acyl_chain_length}, need at least 10 carbons)"
    
    # Next, check the sphingoid “tail”. In our sphingoid match, the last atom (index 3)
    # is expected to connect to a long aliphatic chain.
    sph_tail_start = selected_sphingo[3]
    sph_tail_neighbors = []
    atom_sph_tail = mol.GetAtomWithIdx(sph_tail_start)
    # Exclude the atoms already in the sphingoid core.
    sph_core_set = set(selected_sphingo)
    for nbr in atom_sph_tail.GetNeighbors():
        if nbr.GetIdx() not in sph_core_set and nbr.GetAtomicNum() == 6:
            sph_tail_neighbors.append(nbr.GetIdx())
    if not sph_tail_neighbors:
        return False, "No sphingoid tail (aliphatic chain attached to sphingoid core) found"
    sph_tail_length = 0
    for nbr_idx in sph_tail_neighbors:
        chain_len = longest_carbon_chain(mol, nbr_idx, exclude={sph_tail_start})
        if chain_len > sph_tail_length:
            sph_tail_length = chain_len
    # Set a threshold for the sphingoid tail (e.g. > 6 carbons).
    if sph_tail_length < 6:
        return False, f"Sphingoid tail too short (found chain length {sph_tail_length}, need at least 6 carbons)"
        
    # (Optional) Check overall molecular weight as a rough validation; many ceramides are >500 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for ceramide structure"
        
    return True, ("Found amide bond connecting a long acyl chain (length "
                  f"{acyl_chain_length}) to a phytosphingosine-like backbone with tail length "
                  f"{sph_tail_length} and molecular weight {mol_wt:.1f} Da")
                    
# Example usage:
if __name__ == "__main__":
    # Test the function with one of the provided examples.
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
    valid, reason = is_N_acylphytosphingosine(test_smiles)
    print("Is N-acylphytosphingosine?", valid)
    print("Reason:", reason)