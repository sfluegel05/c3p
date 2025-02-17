"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA

A branched-chain fatty acyl-CoA is defined as a fatty acyl-CoA that results from the condensation of the 
thiol group of coenzyme A with the carboxy group of a branched-chain fatty acid.
This program follows these steps:
  1. Reject SMILES that contain deprotonated groups (e.g. “[O-]”).
  2. Check for a CoA moiety (using an approximate SMARTS match).
  3. Look for thioester groups (SMARTS: "[#6](=O)S") that should be present in acyl-CoAs.
  4. For each thioester:
       a. Identify its carbonyl carbon.
       b. Gather all carbon neighbors of this carbon except those belonging to the oxygen or the sulfur.
          (If more than one neighbor is found, the acyl fragment already splits at the linkage.)
       c. Using DFS that only traverses carbon–carbon bonds (and excluding the carbonyl carbon),
          collect the entire acyl fragment. Also require that the fragment reaches a minimum size and 
          does not include atoms in a ring.
       d. Build the induced subgraph (only carbon–carbon links within the fragment) and check each 
          atom’s degree. In a linear chain, internal carbons should have exactly two neighbors and help 
          detect any branch when an atom has >2.
  5. If any candidate acyl fragment is found to be branched we return True with an explanation.
  
If no candidate thioester yields a branched acyl chain, the function returns False.
"""

from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        Tuple[bool, str]: A boolean indicating if the molecule is a branched-chain fatty acyl-CoA 
                          and a string with the reason for classification.
    """
    # Reject molecules with explicit deprotonation markers.
    if "[O-]" in smiles:
        return False, "Molecule appears deprotonated ([O-] detected); expected neutral CoA form."
    
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check for a CoA motif.
    # (This is an approximate SMARTS for a CoA fragment.)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found."
    
    # Identify thioester groups: a carbonyl carbon bonded to an oxygen (=O) and sulfur.
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found."

    # Helper function: Perform a depth-first search (DFS) starting at a given root (acyl chain starting carbon)
    # and traverse only through carbon atoms. We exclude the carbonyl carbon to avoid backtracking.
    def get_acyl_fragment(root_idx, carbonyl_idx):
        visited = set()
        frontier = [root_idx]
        while frontier:
            curr = frontier.pop()
            if curr in visited:
                continue
            visited.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            # Traverse only neighbors that are carbons.
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == carbonyl_idx:
                    continue  # do not go back into the carbonyl
                if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                    frontier.append(nbr_idx)
        return visited

    # Helper function: Checks for branching in a connected set of carbon atoms.
    # In an unbranched (linear) alkyl chain each carbon (within the induced subgraph) should have 1 (terminal) or 2 (internal) neighbors.
    def is_fragment_branched(fragment_indices):
        for idx in fragment_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Count neighbors that are carbons and lie within the fragment.
            carbon_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() 
                                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in fragment_indices]
            if len(carbon_neighbors) > 2:
                return True
        return False

    # Minimum acyl fragment length (number of carbon atoms) to consider as a fatty acyl chain.
    MIN_FRAGMENT_LENGTH = 3  # adjust as needed

    # For each thioester match, try to extract the fatty acyl fragment and check for branching.
    for match in thioester_matches:
        # match according to "[#6](=O)S": match[0] is the carbonyl carbon, match[1] is the oxygen (from =O), match[2] is sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify candidate acyl root(s): These are carbon neighbors of the carbonyl (excluding the oxygen and sulfur).
        candidate_acyl_roots = []
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip if this neighbor is the specified oxygen or the sulfur.
            if nbr_idx in (match[1], sulfur_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                candidate_acyl_roots.append(nbr_idx)
        
        if not candidate_acyl_roots:
            continue  # cannot extract acyl chain from this thioester candidate
        
        # If more than one candidate acyl root is found, this indicates immediate branching.
        if len(candidate_acyl_roots) > 1:
            explanation = ("Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety. "
                           "Multiple carbon neighbors at the acyl linkage indicate branching.")
            return True, explanation
        
        # Otherwise, take the single candidate.
        acyl_root = candidate_acyl_roots[0]
        # Extract the acyl fragment using DFS.
        fragment = get_acyl_fragment(acyl_root, carbonyl_idx)
        # Require that the fragment is of sufficient length.
        if len(fragment) < MIN_FRAGMENT_LENGTH:
            continue
        # Ensure the fragment is acyclic.
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in fragment):
            continue

        # At this point we now inspect the connectivity of the fragment.
        # In a linear chain, each carbon should have degree 1 (terminal) or 2 (internal) within the fragment.
        if is_fragment_branched(fragment):
            explanation = ("Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety. "
                           "Branching detected in the acyl fragment.")
            return True, explanation
    
    # If none of the thioester candidates showed a branched acyl chain:
    return False, "No branched acyl chain detected in any thioester candidate."


# Example usage:
if __name__ == "__main__":
    # Test with 2-methylbutanoyl-CoA (should be branched)
    example_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(example_smiles)
    print(result, reason)