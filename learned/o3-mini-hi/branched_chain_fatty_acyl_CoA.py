"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA

A branched-chain fatty acyl-CoA is defined as a fatty acyl-CoA where the fatty acyl moiety
contains at least one branch (i.e. a substituent off the "main" chain). We look for:
  (i) a thioester linkage – indicated by a carbonyl bonded to sulfur ([#6](=O)S),
  (ii) a CoA substructure (using an approximate SMARTS match), and
  (iii) branching in the fatty acyl fragment.
  
It is important that we only traverse carbon atoms to extract the acyl fragment and that we
reject molecules in non‐standard protonation states (as indicated by "[O-]" tokens in the SMILES).
Branching is determined by looking at the induced subgraph of the acyl chain – in a linear chain,
each carbon (except terminals) has exactly two carbon neighbors. An atom with more than two neighbors
(in the subgraph) is evidence of branching.
"""

from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Steps:
      1. Reject SMILES that show deprotonated groups (e.g. “[O-]”),
         to avoid misclassification due to non-standard protonation.
      2. Check for a CoA substructure (using an approximate SMARTS).
      3. Identify a thioester group via the SMARTS "[#6](=O)S".
      4. For each thioester, extract the fatty acyl fragment:
            a. Identify the acyl root: the carbon (atomic number 6) bonded to the carbonyl
               (while ignoring bonds to the oxygen and sulfur in the thioester).
            b. Traverse only to carbon neighbors (via a DFS) to collect the candidate acyl fragment.
            c. Reject the candidate if any part of the fragment is in a ring.
      5. In the induced subgraph (only carbons) of the acyl fragment, check every atom’s degree.
         In a linear chain, internal carbons have degree 2 and terminal carbons degree 1.
         If any atom has a degree >2, a branch is present.
      6. If any thioester candidate yields a branched fatty acyl fragment, return True with details;
         otherwise, return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a branched-chain fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Reject molecules with explicit deprotonation markers to avoid non-standard forms.
    if "[O-]" in smiles:
        return False, "Molecule appears deprotonated ([O-] present); expected neutral CoA form."
    
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a CoA motif. This SMARTS is an approximate match to a fragment of Coenzyme A.
    # (This pattern may be refined further if needed.)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"
    
    # Identify a thioester group: a carbonyl carbon connected to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[#6](=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group (acyl-CoA linkage) found."
    
    # Helper to extract the acyl fragment (set of carbon atom indices) from a starting acyl root.
    # We do a DFS that only follows bonds to carbons.
    def get_acyl_fragment(root_idx, carbonyl_idx):
        visited = set()
        frontier = [root_idx]
        while frontier:
            curr = frontier.pop()
            if curr in visited:
                continue
            visited.add(curr)
            atom = mol.GetAtomWithIdx(curr)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Do not go back to the carbonyl carbon
                if nbr_idx == carbonyl_idx:
                    continue
                # Only traverse carbon atoms
                if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                    frontier.append(nbr_idx)
        return visited

    # Helper to check for branching in the acyl fragment.
    # In the induced subgraph, count the number of neighbors (within the fragment) for each atom.
    def is_fragment_branched(fragment_indices):
        for idx in fragment_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Count neighbors that are also in the fragment.
            nbrs_in_fragment = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in fragment_indices]
            # In a linear chain, an internal carbon has exactly 2 carbon neighbors whereas a terminal has 1.
            # Any count >2 indicates branching.
            if len(nbrs_in_fragment) > 2:
                return True
        return False
    
    # For each thioester match, check if its acyl fragment is branched.
    for match in thioester_matches:
        # In the SMARTS "[#6](=O)S", match[0] is the carbonyl carbon,
        # match[1] is the oxygen, and match[2] is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Find the acyl root from the carbonyl: choose a carbon neighbor that is not the oxygen or sulfur.
        acyl_root = None
        for nbr in carbonyl_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in (match[1], sulfur_idx):
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_root = nbr_idx
                break
        if acyl_root is None:
            continue  # try next thioester match
        
        # Extract the acyl fragment from the acyl root (excluding the carbonyl carbon)
        fragment = get_acyl_fragment(acyl_root, carbonyl_idx)
        if len(fragment) < 2:
            continue  # too short to judge
        
        # Ensure that every atom in the fragment is acyclic (i.e. not in any ring)
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in fragment):
            continue  # fragment not suitable
        
        # Check for branching in the fragment:
        # Using the induced subgraph method: if any atom’s number of carbon neighbors (within fragment) >2,
        # then we classify as branched.
        if is_fragment_branched(fragment):
            explanation = "Contains a thioester-linked fatty acyl chain with a branch and a CoA moiety. " \
                          "Branching detected in the acyl fragment."
            return True, explanation
        
        # If no branching was found within the fragment, we do not count this thioester candidate.
    return False, "No branched acyl chain detected in any thioester candidate"


# Example usage:
if __name__ == "__main__":
    # Example: 2-methylbutanoyl-CoA should be classified as branched.
    example_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_branched_chain_fatty_acyl_CoA(example_smiles)
    print(result, reason)