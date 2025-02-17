"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: Long-chain fatty acyl-CoA
Definition:
    A fatty acyl-CoA that results from the condensation of the thiol group of 
    coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA.

    It first checks for a CoA moiety, then looks for a thioester bond 
    such that the S belongs to the CoA fragment. Finally, it attempts to 
    extract the fatty acyl chain (starting at the carbonyl carbon) and counts 
    how many carbon atoms are found. This count (including the carbonyl carbon)
    should be between 13 and 22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Explanation of the classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First attempt to identify the CoA moiety.
    # A typical fragment seen in many fatty acyl-CoAs is:
    #  SCCNC(=O)CCNC(=O)
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_query = Chem.MolFromSmarts(coa_smarts)
    coa_matches = mol.GetSubstructMatches(coa_query)
    if not coa_matches:
        return False, "Coenzyme A moiety not found"
    # Combine all atom indices from the CoA fragment(s)
    coa_idx_set = set()
    for match in coa_matches:
        for idx in match:
            coa_idx_set.add(idx)
            
    # Next, look for a thioester bond.
    # The thioester is characterized by a carbonyl carbon connected to an S atom:
    thioester_smarts = "[#6](=O)[S]"  # simple pattern for a thioester bond
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester_query)
    if not ts_matches:
        return False, "No thioester bond found"

    # We now search for a thioester bond whose S atom is part of the CoA moiety.
    valid_match = None
    for match in ts_matches:
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        if sulfur_idx in coa_idx_set:
            valid_match = (carbonyl_idx, sulfur_idx)
            break
    if valid_match is None:
        return False, "No thioester bond connecting a fatty acid to a CoA moiety found"
    
    carbonyl_idx, sulfur_idx = valid_match
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the neighbor attached to the carbonyl that is the start of the acyl chain.
    # Do not follow the double-bonded oxygen and do not go to the sulfur (which is in CoA).
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # Skip atoms that are not carbon or that are the sulfur linked to CoA.
        if nbr_idx == sulfur_idx:
            continue
        # Also skip oxygen (the carbonyl O)
        if nbr.GetAtomicNum() != 6:
            continue
        acyl_start = nbr_idx
        break
    if acyl_start is None:
        return False, "Fatty acyl chain not found (no carbon neighbor to the carbonyl group)"
    
    # We want to count the number of carbon atoms in the fatty acyl chain,
    # including the carbonyl carbon.
    # We perform a depth-first search (DFS) that only includes carbon atoms
    # and that does not cross into the CoA moiety (which we already have indices for).
    visited = set()
    stack = [carbonyl_idx]
    while stack:
        curr_idx = stack.pop()
        if curr_idx in visited:
            continue
        # Only count carbons
        atom = mol.GetAtomWithIdx(curr_idx)
        if atom.GetAtomicNum() != 6:
            continue
        # Do not enter known CoA atoms (except we already accepted the carbonyl)
        if curr_idx in coa_idx_set and curr_idx != carbonyl_idx:
            continue
        visited.add(curr_idx)
        # Add neighbors that are carbons and are not in CoA
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Allow traversal if neighbor is carbon and not already visited.
            if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                # Also if the neighbor is part of CoA (and not the starting carbonyl), skip it.
                if nbr_idx in coa_idx_set:
                    continue
                stack.append(nbr_idx)
    acyl_chain_length = len(visited)  # count of carbon atoms in the acyl chain
    if acyl_chain_length < 13:
        return False, f"Fatty acyl chain too short (found {acyl_chain_length} carbons; need at least 13)"
    if acyl_chain_length > 22:
        return False, f"Fatty acyl chain too long (found {acyl_chain_length} carbons; need at most 22)"

    # Additional check: Typically fatty acyl-CoAs have high molecular weight (> ~700 Da)
    # due to the CoA part. (This is an optional check.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for fatty acyl-CoA"

    return True, f"Found a valid fatty acyl-CoA with acyl chain length {acyl_chain_length} carbons"

# Sample test (the user provided many SMILES, here just one example is shown)
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    is_valid, reason = is_long_chain_fatty_acyl_CoA(test_smiles)
    print("Is long-chain fatty acyl-CoA:", is_valid)
    print("Reason:", reason)