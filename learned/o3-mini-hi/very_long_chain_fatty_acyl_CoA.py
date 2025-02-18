"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: Very long‐chain fatty acyl‐CoA
Definition: A fatty acyl‐CoA in which the fatty acyl group (the acyl chain attached via a thioester bond)
has a chain length greater than C22 (i.e. at least 23 carbon atoms in the acyl fragment).
Our improved strategy:
  1. Check for a thioester using a SMARTS pattern (C(=O)[S]).
  2. Check for a CoA moiety by looking for an adenine skeleton.
  3. Mark the carbonyl carbon of the thioester and break the bond between it and the sulfur.
  4. Among the resulting fragments, select the one that contains the marked carbon but does NOT contain adenine.
  5. In that fragment, compute the length of the longest contiguous chain of carbons 
     (ignoring dummy atoms and non-carbon atoms).
  6. If the longest “carbon chain” is at least 23, then it qualifies as very long-chain fatty acyl-CoA.
"""

from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    
    Criteria:
       1. Must contain a thioester group (pattern "C(=O)[S]").
       2. Must contain a CoA moiety (detected via adenine substructure patterns).
       3. When the thioester bond (between the carbonyl carbon and sulfur) is broken,
          and the fragment that does not include adenine is examined, the longest contiguous
          carbon chain within that fragment must have at least 23 carbon atoms.
          
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        (bool, str): A tuple of (True, explanation) if the molecule meets the criteria,
                     (False, explanation) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Look for the thioester group: a carbonyl carbon bonded to a sulfur: C(=O)[S]
    thioester_smarts = "C(=O)[S]"
    thioester_pat = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pat is None:
        return False, "Error in creating thioester SMARTS pattern"
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "No thioester group found; not an acyl-CoA"
    # take the first match; indices: carbonyl carbon, then sulfur.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]
    
    # STEP 2: Check for the CoA moiety by looking for adenine substructures.
    adenine_smarts1 = "n1cnc2c(N)ncnc12"
    adenine_smarts2 = "n1cnc2ncnc12"
    adenine_pat1 = Chem.MolFromSmarts(adenine_smarts1)
    adenine_pat2 = Chem.MolFromSmarts(adenine_smarts2)
    if adenine_pat1 is None or adenine_pat2 is None:
        return False, "Error creating adenine SMARTS patterns"
    has_adenine = mol.HasSubstructMatch(adenine_pat1) or mol.HasSubstructMatch(adenine_pat2)
    if not has_adenine:
        return False, "No CoA moiety detected (adenine fragment missing)"
    
    # STEP 3: Mark the carbonyl atom so we can trace its fragment.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    carbonyl_atom.SetProp("is_acyl", "1")
    
    # Find the bond between the carbonyl carbon and the sulfur.
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not found between carbonyl carbon and sulfur"
    bond_idx = bond.GetIdx()
    
    # Fragment the molecule by breaking the thioester bond.
    try:
        frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    except Exception as e:
        return False, f"Error fragmenting molecule: {str(e)}"
    
    # Get fragments as separate molecules.
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments obtained after bond breaking"
    
    # STEP 4: Choose the fragment that contains the marked carbon AND (preferably) does NOT contain adenine.
    acyl_frag = None
    for frag in frags:
        # Check if this fragment contains the marked carbon.
        if any(atom.HasProp("is_acyl") and atom.GetProp("is_acyl")=="1" for atom in frag.GetAtoms()):
            # If this fragment also contains an adenine substructure, it might be part CoA.
            if frag.HasSubstructMatch(adenine_pat1) or frag.HasSubstructMatch(adenine_pat2):
                continue
            acyl_frag = frag
            break
    # Fallback: if none without adenine, choose the fragment that has the mark.
    if acyl_frag is None:
        for frag in frags:
            if any(atom.HasProp("is_acyl") and atom.GetProp("is_acyl")=="1" for atom in frag.GetAtoms()):
                acyl_frag = frag
                break
    if acyl_frag is None:
        return False, "Could not isolate fatty acyl fragment from the thioester bond"
        
    # STEP 5: Define a helper (nested) function to compute the longest contiguous chain of carbon atoms.
    def longest_carbon_chain_length(mol_fragment):
        # Build a graph on the atoms that are carbons.
        # We use a simple DFS for each carbon atom to get the maximum chain length.
        carbon_indices = [atom.GetIdx() for atom in mol_fragment.GetAtoms() if atom.GetAtomicNum() == 6]
        # Build an adjacency dictionary: only count bonds between carbons.
        adj = {idx: [] for idx in carbon_indices}
        for atom in mol_fragment.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            a_idx = atom.GetIdx()
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    adj[a_idx].append(nbr.GetIdx())
                    
        # DFS to find longest simple path in an undirected graph.
        visited_global = set()
        max_length = 0
        
        def dfs(curr, visited):
            local_max = 1
            for nxt in adj.get(curr, []):
                if nxt in visited:
                    continue
                path_length = 1 + dfs(nxt, visited | {nxt})
                if path_length > local_max:
                    local_max = path_length
            return local_max
        
        for start in carbon_indices:
            chain_len = dfs(start, {start})
            if chain_len > max_length:
                max_length = chain_len
        return max_length
    
    # Compute the longest contiguous carbon chain from the acyl fragment.
    lc = longest_carbon_chain_length(acyl_frag)
    if lc < 23:
        return False, f"Longest contiguous carbon chain is {lc} carbons; must be at least 23"
    
    return True, f"Longest contiguous carbon chain in the acyl fragment is {lc} carbons, which is >=23"

# (Optional) Example usage:
# test_smiles = ("CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)"
#                "[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)"
#                "[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
# result, reason = is_very_long_chain_fatty_acyl_CoA(test_smiles)
# print(result, reason)