"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: Sterol Ester
Definition: A steroid ester obtained by formal condensation of the carboxy group 
of any carboxylic acid with the 3-hydroxy group of a sterol.
This implementation:
  1. Checks for at least one ester functional group (SMARTS "C(=O)O").
  2. Identifies rings that are mostly carbons (≥80% C atoms) and then builds fused ring clusters,
     merging rings that share at least 2 atoms.
  3. Keeps only those clusters which have at least 4 rings (with at least one 5-membered ring).
  4. For each ester found, it inspects the ester oxygen (the one single‐bonded to the carbonyl) 
     and its near neighbors (one extra bond away) to see if they belong to any steroid cluster.
If so, the molecule is considered a sterol ester.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.

    A sterol ester is defined as a steroid ester formed by a carboxylic acid with
    the 3-hydroxy group of a sterol. We expect an ester group (C(=O)O) 
    where the alcohol oxygen is connected (directly or via one extra bond)
    to a fused tetracyclic (typically three six-membered and one five-membered ring)
    nucleus that is largely carbon. 
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if classified as a sterol ester, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Look for an ester group: Use SMARTS "C(=O)O"
    ester_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group (C(=O)O) found"
    
    # Step 2. Identify rings that may be part of a steroid nucleus.
    # We consider a ring if at least 80% of its atoms are carbons (for small rings, require at least 4 carbons if 5-membered).
    ring_info = mol.GetRingInfo()
    raw_rings = ring_info.AtomRings()
    # Filter rings for “mostly carbon” rings:
    valid_rings = []
    for ring in raw_rings:
        atom_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        nC = sum(1 for s in atom_symbols if s == "C")
        if len(ring) >= 5:
            if nC/len(ring) >= 0.8:
                valid_rings.append(set(ring))
        else:
            # For very small rings, require at least 4 carbons.
            if nC >= 4:
                valid_rings.append(set(ring))
    if not valid_rings:
        return False, "No rings with steroid-like (largely carbon) character detected"
    
    # Build fused ring clusters.
    # Two rings are fused if they share at least 2 atoms.
    clusters = []
    for ring in valid_rings:
        merged = False
        for cluster in clusters:
            if len(cluster & ring) >= 2:
                cluster |= ring
                merged = True
                break
        if not merged:
            clusters.append(set(ring))
    
    # Now count how many distinct valid rings (from valid_rings) are present in each cluster.
    steroid_clusters = []
    for cluster in clusters:
        count = 0
        has_five = False
        for ring in valid_rings:
            # if at least 2 atoms of this ring appear in the cluster, assume it is represented.
            if len(ring & cluster) >= 2:
                count += 1
                if len(ring) == 5:
                    has_five = True
        # A typical steroid nucleus has 4 fused rings (with one being five-membered)
        if count >= 4 and has_five:
            steroid_clusters.append(cluster)
    if not steroid_clusters:
        return False, ("Steroid nucleus not detected – expected a fused system of at least 4 rings "
                       "(with at least one 5-membered ring and largely carbon in nature)")
    
    # Step 3. For each ester match, locate its candidate ester oxygen (the one bound by a single bond to the C=O).
    sterol_ester_found = False
    for match in ester_matches:
        # Identify the carbonyl carbon and candidate oxygens in the matching substructure.
        candidate_carbon = None
        candidate_oxygens = []
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                candidate_carbon = idx
            elif atom.GetAtomicNum() == 8:
                candidate_oxygens.append(idx)
        if candidate_carbon is None or not candidate_oxygens:
            continue
        
        # For each candidate oxygen, check that the bond to the carbonyl carbon is a single bond.
        for oxy_idx in candidate_oxygens:
            bond = mol.GetBondBetweenAtoms(candidate_carbon, oxy_idx)
            if bond is None or bond.GetBondTypeAsDouble() != 1.0:
                continue
            # Now, determine the neighborhood of this oxygen.
            # We want to see if it (or a neighbor up to one extra bond away) is a member of a steroid cluster.
            candidate_set = {oxy_idx}
            # First-degree neighbors (except the carbonyl carbon)
            first_neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(oxy_idx).GetNeighbors() if nbr.GetIdx() != candidate_carbon]
            candidate_set.update(first_neighbors)
            # Include second-degree neighbors (neighbors of the first neighbors)
            second_neighbors = set()
            for n_idx in first_neighbors:
                for nn in mol.GetAtomWithIdx(n_idx).GetNeighbors():
                    if nn.GetIdx() not in candidate_set:
                        second_neighbors.add(nn.GetIdx())
            candidate_set.update(second_neighbors)
            
            # Check overlap with any steroid cluster.
            for cluster in steroid_clusters:
                if candidate_set & cluster:
                    sterol_ester_found = True
                    break
            if sterol_ester_found:
                break
        if sterol_ester_found:
            break

    if not sterol_ester_found:
        return False, ("Ester group found but the candidate ester oxygen (or its near neighbors) is not "
                       "connected to a fused steroid-like nucleus (expected for sterol-derived alcohol)")
    
    return True, ("Contains at least one ester group where the alcohol oxygen (directly or via a short bond) "
                  "is connected to a fused tetracyclic steroid nucleus (with at least one 5-membered ring), "
                  "consistent with a sterol ester")

# Example usage:
if __name__ == "__main__":
    # Try an example SMILES from the provided list.
    smiles_example = "O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O"
    result, reason = is_sterol_ester(smiles_example)
    print("Result:", result)
    print("Reason:", reason)