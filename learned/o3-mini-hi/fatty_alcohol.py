"""
Classifies: CHEBI:24026 fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: fatty alcohol
Definition: "An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
             Fatty alcohols may be saturated or unsaturated and may be branched or unbranched."

This implementation:
  1. Parses the molecule using RDKit.
  2. Searches for alcohol groups defined as a saturated (or stereochemically annotated) carbon attached to –OH.
  3. For each candidate –OH we require that the alcohol carbon is not in a ring and its non‐OH neighbors are carbon.
  4. We then build a graph of acyclic, non‐aromatic carbon atoms (including sp2 carbons from unsaturated chains).
  5. We “grow” the longest simple (non‐repeating) chain (using DFS) starting from the alcohol carbon.
  6. To avoid “contaminated” chains we also require that each carbon in the candidate chain does not have an extra substituent
     that is a heteroatom.
  7. If one candidate –OH gives a longest chain of at least MIN_CHAIN_LENGTH (here 7) carbons, we classify the molecule as fatty alcohol.
  
Note:
  We no longer reject molecules with other oxygenated (or carboxylate) groups.
"""

from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule belongs to the fatty alcohol class.
    
    The molecule must:
      - Parse correctly.
      - Contain at least one alcohol group defined as an sp3 (or stereochemically annotated) carbon bound to an -OH.
      - The candidate alcohol carbon must be acyclic (not in a ring) and, apart from the hydroxyl oxygen,
        be bound only to carbons.
      - From that alcohol carbon, one must be able to follow a contiguous chain of non‐aromatic, non‐ring carbon atoms
        (allowing unsaturation) that is “pure” (i.e. none of the carbons in that chain have extra connections to heteroatoms)
        and has at least MIN_CHAIN_LENGTH atoms.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple; True if classified as a fatty alcohol, False otherwise,
                   plus a human‐readable reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for an alcohol group: a saturated (or annotated) carbon attached to a hydroxyl.[1]
    # (This will catch [CX4] or stereochemically specified [C@H] or [C@@H] when bonded to [OX2H])
    alcohol_smarts = "[C;!R][OX2H]"
    alcohol_query = Chem.MolFromSmarts(alcohol_smarts)
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    if not alcohol_matches:
        return False, "No proper alcohol (-OH) group found on a non‐ring carbon"
    
    # Build a graph (dictionary) of carbons that are aliphatic (non-aromatic and not in a ring)
    # We include sp3 and non-aromatic sp2 carbons, so that unsaturated fatty chains count.
    carbon_nodes = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and (not atom.IsInRing()):
            carbon_nodes.add(atom.GetIdx())
    # Build a connectivity graph: only between carbons in carbon_nodes.
    carbon_graph = {idx: [] for idx in carbon_nodes}
    for idx in carbon_nodes:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in carbon_nodes:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # DFS function that returns (length, path) for the longest simple path starting from 'current'.
    def dfs_longest(current, visited):
        best_path = [current]
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                path = dfs_longest(nbr, visited | {nbr})
                if len(path) + 1 > len(best_path):
                    best_path = [current] + path
        return best_path

    # Helper function to check “purity” of a candidate chain:
    # For each carbon in the chain, besides its chain neighbors, it should not have attached heteroatoms.
    def chain_is_pure(chain):
        for cid in chain:
            atom = mol.GetAtomWithIdx(cid)
            chain_neighbors = set()
            # In the chain, consider next and previous (if any)
            # (We check only the bonds along the chain.)
            for other in chain:
                if other == cid:
                    continue
                # If atoms are directly bonded, add them.
                if mol.GetAtomWithIdx(cid).HasBondBetween(mol.GetAtomWithIdx(other)):
                    chain_neighbors.add(other)
            # Now examine all neighbors of this atom.
            for nbr in atom.GetNeighbors():
                # If neighbor is not in the chain and is not hydrogen,
                # then if it is not carbon then the chain is “decorated” with a heteroatom.
                if nbr.GetAtomicNum() != 1 and nbr.GetIdx() not in chain_neighbors:
                    if nbr.GetAtomicNum() != 6:
                        return False
        return True


    # Use a minimum chain length requirement.
    MIN_CHAIN_LENGTH = 7  # counting the alcohol carbon itself
    
    best_chain_length = 0
    candidate_found = False
    candidate_reason = ""
    
    # Loop over candidate alcohol groups.
    for match in alcohol_matches:
        # match is a tuple: first element is the candidate alcohol carbon and the second the hydroxyl oxygen.
        alc_c_idx = match[0]
        alc_o_idx = match[1]
        alc_atom = mol.GetAtomWithIdx(alc_c_idx)
        # Skip if the alcohol carbon is in a ring (we want acyclic chain)
        if alc_atom.IsInRing():
            continue
        
        # Check that aside from the -OH, the alcohol carbon's neighbors are all carbons.
        valid = True
        for nbr in alc_atom.GetNeighbors():
            if nbr.GetIdx() == alc_o_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                valid = False
                break
        if not valid:
            continue
        
        # Check that the candidate alcohol carbon is in our aliphatic carbon graph.
        if alc_c_idx not in carbon_graph:
            continue
        
        # Compute the longest chain (simple path) from the alcohol carbon.
        longest_path = dfs_longest(alc_c_idx, {alc_c_idx})
        chain_length = len(longest_path)
        if chain_length > best_chain_length:
            best_chain_length = chain_length
        
        # Check chain “purity”: along the chain, no carbon should have extra hetero-atom substituents.
        if not chain_is_pure(longest_path):
            # Skip this candidate if the chain is "decorated" by non-carbon atoms.
            continue
        
        if chain_length >= MIN_CHAIN_LENGTH:
            candidate_found = True
            candidate_reason = (f"Molecule contains an alcohol group attached to a "
                                f"{chain_length}-carbon acyclic aliphatic chain (chain path: {longest_path}).")
            break  # accept the first acceptable candidate

    if candidate_found:
        return True, candidate_reason
    else:
        return False, (f"Longest acyclic aliphatic chain from any candidate alcohol group is only "
                       f"{best_chain_length} carbons (need at least {MIN_CHAIN_LENGTH}).")
        
# Example usage (for testing):
if __name__ == "__main__":
    # A true positive example: octan-2-ol (SMILES: CCCCCCC(C)O)
    test_smiles = [
        "CCCCCCC(C)O",                       # octan-2-ol (should be true)
        "CCCC[C@H](C)O",                      # (2S)-2-heptanol (borderline; our minimum chain length is 7)
        "OC(CCCCCCCCCCCCCCC)CC(=O)C1=CC=CC=C1",# 3-Hydroxy-1-phenyl-1-octadecanone
        "OC1=CC(O)=CC(=C1)CCCCCCCCCCCCCCC[C@@H](O)C", # Hansfordiol A
        "O=C(O[C@@H](CCCCCCCCCCCCC[C@H](O)C)[C@H](C1=C)C(=O)O"  # Allo-murolic acid (has extra acid groups)
    ]
    for smi in test_smiles:
        result, reason = is_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nClassified as fatty alcohol? {result}\nReason: {reason}\n{'-'*60}")