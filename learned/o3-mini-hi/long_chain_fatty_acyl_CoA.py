"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: Long-chain fatty acyl-CoA
Definition:
    A fatty acyl-CoA that results from the condensation of the thiol group of 
    coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.
    
This program uses an improved approach:
  • It requires the presence of a thioester bond (C(=O)–S) 
    connecting the fatty acyl chain to a CoA moiety.
  • Instead of hardcoding the full CoA fragment, it looks for an adenine 
    substructure (a key component of CoA).
  • It then verifies that the sulfur atom of the thioester can reach 
    adenine through the molecular graph.
  • Finally, it extracts the acyl chain (by DFS from the carbonyl carbon, 
    avoiding the sulfur side) and counts its carbon atoms (should be within 13–22).
    
If any of these conditions fail, a descriptive reason is returned.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA.
    
    Steps:
      1. Parse the molecule.
      2. Look for an adenine substructure (indicative of a CoA moiety).
      3. Look for a thioester bond: a carbonyl carbon ([#6](=O)) bonded to a sulfur ([S]).
      4. Verify that the sulfur atom connects (via non-carbonyl paths) to the adenine set.
      5. From the carbonyl carbon, traverse only carbon atoms (avoiding the sulfur side)
         to count how many atoms belong to the fatty acyl chain.
      6. Ensure that the chain length (including the carbonyl carbon) is between 13 and 22.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Look for adenine substructure as a proxy for the CoA moiety.
    # A common SMARTS for adenine (allowing for substituents) is:
    adenine_smarts = "n1cnc2c1ncnc2"  
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    adenine_matches = mol.GetSubstructMatches(adenine_query)
    if not adenine_matches:
        return False, "Coenzyme A moiety not found (adenine fragment missing)"
    # Union all atoms in adenine matches into a set
    adenine_idx_set = set()
    for match in adenine_matches:
        adenine_idx_set.update(match)
    
    # Step 2: Find thioester bonds.
    # The pattern [#6](=O)[S] finds a carbonyl carbon bonded to sulfur.
    thioester_smarts = "[#6](=O)[S]"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester_query)
    if not ts_matches:
        return False, "No thioester bond found"
    
    valid_match = None
    # For each thioester bond, check if the sulfur connects to adenine.
    for match in ts_matches:
        carbonyl_idx, sulfur_idx = match  # match: (carbonyl, sulfur)
        # Perform a DFS starting from the sulfur atom, avoiding the carbonyl atom.
        stack = [sulfur_idx]
        visited = set()
        found_adenine = False
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            if current in adenine_idx_set:
                found_adenine = True
                break
            # Add all neighbors except the carbonyl atom
            for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == carbonyl_idx:
                    continue
                if nbr_idx not in visited:
                    stack.append(nbr_idx)
        if found_adenine:
            valid_match = (carbonyl_idx, sulfur_idx)
            break
    if valid_match is None:
        return False, "No thioester bond connecting a fatty acid to a CoA moiety found"
    carbonyl_idx, sulfur_idx = valid_match
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Step 3: From the carbonyl carbon, determine the acyl chain.
    # We identify the neighbor attached to the carbonyl that is not the sulfur or the carbonyl oxygen.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # Skip the sulfur atom in the thioester bond
        if nbr_idx == sulfur_idx:
            continue
        # Skip non-carbon atoms (e.g. carbonyl oxygen)
        if nbr.GetAtomicNum() != 6:
            continue
        acyl_start = nbr_idx
        break
    if acyl_start is None:
        return False, "Fatty acyl chain not found (no carbon neighbor to the carbonyl group)"
    
    # Step 4: Traverse from the carbonyl atom (only along carbon atoms)
    # using a DFS that does not cross over to the CoA side.
    # We will not cross back through the sulfur of the thioester bond.
    acyl_visited = set()
    stack = [carbonyl_idx]
    while stack:
        current = stack.pop()
        if current in acyl_visited:
            continue
        # Only consider carbon atoms
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() != 6:
            continue
        acyl_visited.add(current)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Avoid going into the thioester bond that leads to CoA
            if nbr_idx == sulfur_idx:
                continue
            # Only traverse carbons
            if nbr.GetAtomicNum() == 6 and nbr_idx not in acyl_visited:
                stack.append(nbr_idx)
    acyl_chain_length = len(acyl_visited)  # includes the carbonyl carbon
    
    # Step 5: Check chain length (should be between 13 and 22 carbons).
    if acyl_chain_length < 13:
        return False, f"Fatty acyl chain too short (found {acyl_chain_length} carbons; need at least 13)"
    if acyl_chain_length > 22:
        return False, f"Fatty acyl chain too long (found {acyl_chain_length} carbons; need at most 22)"
    
    # Optional: Check molecular weight (fatty acyl-CoAs are typically high, > ~700 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is lower than expected for fatty acyl-CoA"
    
    return True, f"Found a valid long-chain fatty acyl-CoA with acyl chain length {acyl_chain_length} carbons"

# Sample test: (You can test with any of the provided SMILES.)
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    is_valid, reason = is_long_chain_fatty_acyl_CoA(test_smiles)
    print("Is long-chain fatty acyl-CoA:", is_valid)
    print("Reason:", reason)