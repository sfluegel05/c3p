"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
#!/usr/bin/env python
"""
Classifies: Long-chain fatty acyl-CoA

Definition:
  A fatty acyl-CoA that results from the condensation of the thiol group of coenzyme A with the 
  carboxy group of any long-chain (C13 to C22) fatty acid.

This implementation:
  • Looks for an adenine fragment as a marker for CoA.
  • Searches for a thioester bond (carbonyl attached to sulfur) using SMARTS.
  • Verifies that the sulfur is connected (via DFS) to the adenine portion.
  • Finds the acyl chain starting at the carbonyl carbon (ignoring the S-side) and then 
    enumerates all simple (acyclic) carbon-only paths (i.e. the “linear” chains) to determine the 
    maximum length.
  • Finally, it requires that the chain (counting the carbonyl carbon) be between 13 and 22 carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a long-chain fatty acyl-CoA, else False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Look for an adenine substructure (as a part of coenzyme A)
    adenine_smarts = "n1cnc2c1ncnc2"
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    adenine_matches = mol.GetSubstructMatches(adenine_query)
    if not adenine_matches:
        return False, "Coenzyme A moiety not found (adenine fragment missing)"
    adenine_idx_set = set()
    for match in adenine_matches:
        adenine_idx_set.update(match)
    
    # Step 2: Locate a thioester bond.
    # We use a strict SMARTS pattern: a carbonyl carbon ([CX3]=O) directly bonded to a sulfur.
    thioester_smarts = "[CX3](=O)[S]"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    ts_matches = mol.GetSubstructMatches(thioester_query)
    if not ts_matches:
        return False, "No thioester bond found"
    
    valid_match = None
    # For each thioester match, retrieve the carbonyl carbon (first atom) and the sulfur (third atom).
    for match in ts_matches:
        if len(match) < 3:
            continue
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        # Step 3: From sulfur, follow bonds (excluding the carbonyl side) to see if adenine is reached.
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
    # Select the neighbor (of the carbonyl) that is a carbon and not on the sulfur side.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() == sulfur_idx:
            continue
        if nbr.GetAtomicNum() == 6:
            acyl_start = nbr.GetIdx()
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found (no carbon neighbor to the carbonyl group)"
    
    # Step 5: Enumerate all simple paths obtained by walking only through carbon atoms.
    # We include the carbonyl carbon as the first atom in our chain.
    # Since fatty acyl chains are linear (acyclic) backbones, we enumerate only simple paths.
    def dfs_paths(current_path):
        """
        Recursively enumerate all simple paths (acyclic) that extend the current carbon chain.
        current_path is a list of atom indices.
        Returns a list of paths (each a list of indices) ending at a terminal carbon 
        (no further carbon neighbors that can extend the path).
        """
        paths = []
        current_atom = mol.GetAtomWithIdx(current_path[-1])
        extended = False
        for nbr in current_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Only follow carbon neighbors that are not already in the path.
            if nbr.GetAtomicNum() == 6 and nbr_idx not in current_path:
                # To avoid “backtracking” to the carbonyl at the start, we don’t allow that either.
                if nbr_idx == carbonyl_idx:
                    continue
                new_path = current_path + [nbr_idx]
                # Recurse from this new path.
                sub_paths = dfs_paths(new_path)
                if sub_paths:
                    paths.extend(sub_paths)
                    extended = True
        if not extended:
            # No further carbon neighbors, so return the current terminal path.
            return [current_path]
        return paths

    # Start the path from [carbonyl, acyl_start]
    initial_path = [carbonyl_idx, acyl_start]
    all_paths = dfs_paths(initial_path)
    if not all_paths:
        return False, "No acyl chain could be traversed"
    # Compute the length (number of carbons in the chain) for every path.
    # The chain length is defined including the carbonyl carbon.
    max_chain_length = max(len(path) for path in all_paths)
    
    # Step 6: Validate chain length between 13 and 22 carbons.
    if max_chain_length < 13:
        return False, f"Fatty acyl chain too short (found {max_chain_length} carbons; need at least 13)"
    if max_chain_length > 22:
        return False, f"Fatty acyl chain too long (found {max_chain_length} carbons; need at most 22)"
    
    # Optionally, check molecular weight (typical fatty acyl-CoAs are >700 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight ({mol_wt:.1f} Da) lower than expected for fatty acyl-CoA"
    
    return True, f"Found a valid long-chain fatty acyl-CoA with acyl chain length {max_chain_length} carbons"

# Sample test (feel free to test with provided SMILES examples)
if __name__ == "__main__":
    # A sample SMILES string for hexadec-2-enoyl-CoA from the given list.
    test_smiles = "[H]C(CCCCCCCCCCCCC)=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    valid, reason = is_long_chain_fatty_acyl_CoA(test_smiles)
    print("Is long-chain fatty acyl-CoA:", valid)
    print("Reason:", reason)