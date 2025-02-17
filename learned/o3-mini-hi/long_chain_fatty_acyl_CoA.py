"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: Long-chain fatty acyl-CoA

Definition:
  A fatty acyl-CoA that results from the condensation of the thiol group of coenzyme A with the 
  carboxy group of any long-chain (C13 to C22) fatty acid.

This approach:
  • Finds an adenine substructure (as a marker for CoA).
  • Locates a thioester bond using “[#6](=O)[S]” (i.e. the carbonyl carbon attached to oxygen then sulfur).
  • Verifies that the sulfur side connects (via a DFS) to an adenine atom.
  • From the carbonyl carbon it picks the carbon neighbor (excluding the sulfur) as the acyl chain start.
  • It then “walks” the chain via carbon–carbon bonds using a recursive search that tracks visited atoms.
  • Finally it checks that the chain length (counting the carbonyl carbon) is between 13 and 22.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Look for an adenine substructure (as a part of CoA).
    adenine_smarts = "n1cnc2c1ncnc2"
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    adenine_matches = mol.GetSubstructMatches(adenine_query)
    if not adenine_matches:
        return False, "Coenzyme A moiety not found (adenine fragment missing)"
    adenine_idx_set = set()
    for match in adenine_matches:
        adenine_idx_set.update(match)
    
    # Step 2: Locate a thioester bond. SMARTS: carbonyl carbon bonded to oxygen then to sulfur.
    # Note: The SMARTS "[#6](=O)[S]" matches a carbonyl (C=O) directly bonded to an S.
    thioester_smarts = "[#6](=O)[S]"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester_query)
    if not ts_matches:
        return False, "No thioester bond found"
    
    valid_match = None
    # For each thioester match, get the carbonyl carbon (first atom) and the sulfur (third atom).
    for match in ts_matches:
        # match is a tuple of indices corresponding to the atoms in the SMARTS
        if len(match) < 3:
            continue
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        # Step 3: Verify that starting from the sulfur (excluding the carbonyl side)
        # we can eventually reach an adenine atom.
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
                if nbr_idx == carbonyl_idx:
                    continue
                if nbr_idx not in visited:
                    stack.append(nbr_idx)
        if found_adenine:
            valid_match = (carbonyl_idx, sulfur_idx)
            break
    if valid_match is None:
        return False, "No thioester bond connecting fatty acid to CoA moiety found"
    
    carbonyl_idx, sulfur_idx = valid_match
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Step 4: Identify the acyl chain start.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # Do not follow the sulfur side.
        if nbr_idx == sulfur_idx:
            continue
        if nbr.GetAtomicNum() == 6:
            acyl_start = nbr_idx
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found (no carbon neighbor to the carbonyl group)"
    
    # Step 5: Traverse the acyl chain along carbon–carbon bonds.
    # To avoid infinite recursion due to cycles, we include a visited set.
    def longest_linear_chain(curr_idx, visited):
        """
        Recursively computes the length of the longest linear chain starting from curr_idx.
        visited is a set of atom indices already traversed to avoid cycles.
        Returns the number of carbons in the chain starting at curr_idx.
        """
        max_length = 1  # count the current carbon
        # Get the current atom.
        curr_atom = mol.GetAtomWithIdx(curr_idx)
        for nbr in curr_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Only follow carbon–carbon bonds that haven't been visited.
            if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                # Create a new visited set for the recursive branch.
                new_visited = set(visited)
                new_visited.add(nbr_idx)
                branch_length = 1 + longest_linear_chain(nbr_idx, new_visited)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    try:
        # Start the chain by including the carbonyl group.
        # The acyl chain length is defined as the carbonyl carbon plus the chain.
        initial_visited = {carbonyl_idx, acyl_start}
        chain_length = 1 + longest_linear_chain(acyl_start, initial_visited)
    except RecursionError:
        return False, "Recursion error during acyl chain traversal"

    # Step 6: Validate chain length between 13 and 22 carbons.
    if chain_length < 13:
        return False, f"Fatty acyl chain too short (found {chain_length} carbons; need at least 13)"
    if chain_length > 22:
        return False, f"Fatty acyl chain too long (found {chain_length} carbons; need at most 22)"
    
    # Optional: Check molecular weight (typical fatty acyl-CoAs are >700 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight ({mol_wt:.1f} Da) lower than expected for fatty acyl-CoA"
    
    return True, f"Found a valid long-chain fatty acyl-CoA with acyl chain length {chain_length} carbons"

# Sample test (feel free to test with provided SMILES examples)
if __name__ == "__main__":
    test_smiles = "CCCCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    valid, reason = is_long_chain_fatty_acyl_CoA(test_smiles)
    print("Is long-chain fatty acyl-CoA:", valid)
    print("Reason:", reason)