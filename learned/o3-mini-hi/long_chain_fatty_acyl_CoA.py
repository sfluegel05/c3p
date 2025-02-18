"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: Long-chain fatty acyl-CoA
Definition:
    A fatty acyl-CoA that results from the condensation of the thiol group of 
    coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.
    
This code:
  • Searches for an adenine substructure (a key component of CoA).
  • Finds a thioester bond via the SMARTS pattern "[#6](=O)[S]". Note that this pattern returns three atoms; 
    we use the first (carbonyl carbon) and third (sulfur) atoms.
  • Verifies that the thioester bond connects to the CoA moiety (by ensuring that starting from the sulfur 
    we can reach an adenine atom).
  • Traverses the carbon chain from the carbonyl carbon (avoiding the sulfur) and counts contiguous carbon atoms,
    which is expected to be between 13 and 22.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA.
    
    Steps:
      1. Parse the molecule.
      2. Look for an adenine substructure (indicative of a CoA moiety).
      3. Look for a thioester bond using the pattern "[#6](=O)[S]". (Unpack the carbonyl carbon and sulfur.)
      4. Verify that a DFS starting from the sulfur (avoiding the carbonyl carbon) reaches an adenine atom.
      5. From the carbonyl carbon, traverse only along carbon atoms (avoiding the sulfur side) to identify the fatty acyl chain.
      6. Check that the chain length (including the carbonyl carbon) is between 13 and 22 carbons.
    
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
    
    # Step 1: Look for adenine substructure to indicate the CoA moiety.
    # SMARTS for adenine (allowing for substituents)
    adenine_smarts = "n1cnc2c1ncnc2"  
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    adenine_matches = mol.GetSubstructMatches(adenine_query)
    if not adenine_matches:
        return False, "Coenzyme A moiety not found (adenine fragment missing)"
    adenine_idx_set = set()
    for match in adenine_matches:
        adenine_idx_set.update(match)
    
    # Step 2: Find thioester bonds.
    # The pattern "[#6](=O)[S]" returns three atoms: carbonyl carbon, oxygen, and sulfur.
    thioester_smarts = "[#6](=O)[S]"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester_query)
    if not ts_matches:
        return False, "No thioester bond found"
    
    valid_match = None
    # For each thioester bond, retrieve the carbonyl carbon and the sulfur (skip the oxygen).
    for match in ts_matches:
        # match will contain 3 atoms: [carbonyl, oxygen, sulfur]
        if len(match) < 3:
            continue  # Safety check, though should always have 3 atoms
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        # Step 3: Check if the sulfur connects (via non-carbonyl paths) to adenine.
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
            for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Avoid going back to the carbonyl carbon.
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
    
    # Step 4: Find the starting atom of the fatty acyl chain.
    # It must be a carbon neighbor of the carbonyl but not the sulfur.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        if nbr_idx == sulfur_idx:  # Skip the thioester bond side that leads to CoA
            continue
        if nbr.GetAtomicNum() != 6:  # Only consider carbon atoms
            continue
        acyl_start = nbr_idx
        break
    if acyl_start is None:
        return False, "Fatty acyl chain not found (no carbon neighbor to the carbonyl group)"
    
    # Step 5: Traverse the acyl chain.
    # We start from the carbonyl carbon and only traverse via carbon atoms.
    acyl_visited = set()
    stack = [carbonyl_idx]
    while stack:
        current = stack.pop()
        if current in acyl_visited:
            continue
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() != 6:
            continue
        acyl_visited.add(current)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Avoid going through the thioester linkage toward CoA.
            if nbr_idx == sulfur_idx:
                continue
            # Only traverse further if the neighbor is carbon.
            if nbr.GetAtomicNum() == 6 and nbr_idx not in acyl_visited:
                stack.append(nbr_idx)
    acyl_chain_length = len(acyl_visited)  # This count includes the carbonyl carbon.
    
    # Step 6: Validate chain length (should be between 13 and 22 carbons).
    if acyl_chain_length < 13:
        return False, f"Fatty acyl chain too short (found {acyl_chain_length} carbons; need at least 13)"
    if acyl_chain_length > 22:
        return False, f"Fatty acyl chain too long (found {acyl_chain_length} carbons; need at most 22)"
    
    # Optional check: molecular weight for fatty acyl-CoAs is typically high (>700 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight ({mol_wt:.1f} Da) lower than expected for fatty acyl-CoA"
    
    return True, f"Found a valid long-chain fatty acyl-CoA with acyl chain length {acyl_chain_length} carbons"

# Sample test (you can test with different provided SMILES)
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    is_valid, reason = is_long_chain_fatty_acyl_CoA(test_smiles)
    print("Is long-chain fatty acyl-CoA:", is_valid)
    print("Reason:", reason)