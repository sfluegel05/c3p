"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: Sterol Ester
Definition: A steroid ester obtained by formal condensation of the carboxy group 
of any carboxylic acid with the 3-hydroxy group of a sterol.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.

    A sterol ester is defined as a steroid ester formed by reaction of any 
    carboxylic acid with the 3-hydroxy group of a sterol. Thus we expect
    an ester group (C(=O)O) where the oxygen from the sterol (the alcohol side)
    is connected (directly or through one extra bond) to a fused tetracyclic 
    ring system (three six-membered rings and one five-membered ring).

    The algorithm first: 
      1. Confirms the presence of at least one ester (using SMARTS C(=O)O).
      2. Identifies fused ring clusters (by merging rings that share at least 2 atoms)
         and then retains clusters that have at least four rings, with one being 5-membered.
      3. For each ester group found, it examines the candidate “ester oxygen” (the oxygen
         attached by a single bond to the carbonyl carbon) and checks whether either it,
         or an atom one bond away from it, is in one such steroid-like fused ring cluster.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a sterol ester, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1. Look for an ester functional group using the pattern "C(=O)O"
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group (C(=O)O) found"
    
    # STEP 2. Identify all rings.
    ring_info = mol.GetRingInfo()
    rings = [set(r) for r in ring_info.AtomRings()]
    if not rings:
        return False, "No rings found in molecule; steroid nucleus expected"
    
    # Build fused ring clusters.
    # Two rings are considered fused if they share at least 2 atoms.
    # We merge rings into clusters by a union-find like approach.
    clusters = []
    for ring in rings:
        added = False
        for clust in clusters:
            # If intersection has 2 or more atoms, merge ring into this cluster.
            if len(clust & ring) >= 2:
                clust |= ring
                added = True
                break
        if not added:
            clusters.append(set(ring))
    # NOTE: The clusters here are sets of atom indices that are in fused rings.
    # However, we do not yet know how many distinct rings each cluster corresponds to.
    # So we will reassemble the rings into clusters: for each cluster, count how many
    # of our original rings (from 'rings') are completely or partly contained.
    steroid_clusters = []
    for clust in clusters:
        count_rings = 0
        has_5_membered = False
        for ring in rings:
            # if at least 50% of the ring's atoms are in cluster, count it
            if len(ring & clust) >= 2:  # using 2 as threshold
                count_rings += 1
                if len(ring) == 5:
                    has_5_membered = True
        # Typical steroid nucleus should have 4 fused rings (3 six-membered and 1 five-membered)
        if count_rings >= 4 and has_5_membered:
            steroid_clusters.append(clust)
    
    if not steroid_clusters:
        return False, ("Steroid nucleus not detected – expected a fused system of at least 4 rings "
                       "(with one 5-membered ring)")
    
    # STEP 3. For each ester match, try to identify the ester oxygen from the alcohol side.
    # The SMARTS "C(=O)O" will match the carbonyl group (C and O attached).
    # The match tuple may include both the carbonyl carbon and the oxygen.
    # For each such match, find the oxygen that is single-bonded to the carbonyl carbon.
    sterol_ester_found = False
    for match in ester_matches:
        # Identify the carbonyl carbon and candidate oxygens.
        candidate_c = None
        candidate_oxygens = []
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                candidate_c = idx
            elif atom.GetAtomicNum() == 8:
                candidate_oxygens.append(idx)
        if candidate_c is None or not candidate_oxygens:
            continue
        
        for oxy_idx in candidate_oxygens:
            bond = mol.GetBondBetweenAtoms(candidate_c, oxy_idx)
            if bond is None:
                continue
            # Accept only if this bond is a single bond.
            if bond.GetBondTypeAsDouble() != 1.0:
                continue
            # candidate ester oxygen from the alcohol side is oxy_idx.
            # Check its direct neighbors (other than the carbonyl carbon).
            candidate_atoms = {oxy_idx}
            neighbor_idxs = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(oxy_idx).GetNeighbors() if nbr.GetIdx() != candidate_c]
            candidate_atoms.update(neighbor_idxs)
            # Also, include second neighbors (i.e. atoms one bond further) to allow a short intervening bond.
            second_neighbors = set()
            for n_idx in neighbor_idxs:
                atom_n = mol.GetAtomWithIdx(n_idx)
                for nn in atom_n.GetNeighbors():
                    if nn.GetIdx() not in candidate_atoms:
                        second_neighbors.add(nn.GetIdx())
            candidate_atoms.update(second_neighbors)
            
            # Check if any candidate atom belongs to one of the steroid clusters.
            for cluster in steroid_clusters:
                if candidate_atoms & cluster:
                    sterol_ester_found = True
                    break
            if sterol_ester_found:
                break
        if sterol_ester_found:
            break
    
    if not sterol_ester_found:
        return False, ("Ester group found but the candidate ester oxygen (or its near neighbors) is not "
                       "connected to a fused steroid-like nucleus (expected for a sterol-derived alcohol)")
    
    return True, ("Contains at least one ester group where the alcohol oxygen (directly or via a short bond) "
                  "is connected to a fused tetracyclic steroid nucleus (with at least one 5-membered ring), "
                  "consistent with a sterol ester")


# Example usage:
if __name__ == "__main__":
    # Try one of the examples provided.
    smiles_example = "O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O"
    result, reason = is_sterol_ester(smiles_example)
    print("Result:", result)
    print("Reason:", reason)