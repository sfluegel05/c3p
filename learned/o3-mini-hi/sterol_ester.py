"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: Sterol Ester
Definition: A steroid ester obtained by formal condensation of the carboxy group 
of any carboxylic acid with the 3-hydroxy group of a sterol.

In this implementation:
  1. We first look for an ester functional group “C(=O)O”
  2. For each ester we identify the “alcohol oxygen” (that is not the carbonyl oxygen) by checking
     that its bond to the carbonyl carbon is single.
  3. For each such candidate oxygen we inspect its (non–carbonyl) neighbor(s) that should be part of
     the sterol. We gather all rings (from mol.GetRingInfo()) that contain this neighbor and that are non‐aromatic.
  4. We “cluster” these rings – that is, we group rings that are fused (we require at least 2 atoms in common).
  5. If any cluster has at least 4 rings and at least one ring is 5‐membered (a typical steroid nucleus), 
     then the ester oxygen is considered to derive from a sterol.
  6. In that case we classify the molecule as a sterol ester.
  
If any step fails, we return a diagnosis as to why.
  
Note: This heuristic is not perfect – many borderline cases exist.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.

    A sterol ester is defined as a steroid ester formed by the condensation of a carboxylic acid 
    and the 3-hydroxy group of a sterol. We expect a C(=O)O ester group where the alcohol oxygen is 
    (either directly or via one extra bond) connected to a fused tetracyclic (or near‐tetracyclic) 
    non‐aromatic ring system containing at least 4 fused rings with at least one 5‐membered ring.

    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if classified as a sterol ester, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Look for an ester group using a simple SMARTS.
    # This pattern matches a carbonyl carbon (C) double-bonded to O and also bonded to an oxygen.
    ester_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group (C(=O)O) found"

    # Helper: given a list of rings (each a set of atom indices), cluster them
    # if they share at least 2 atoms. We use a simple iterative merging scheme.
    def cluster_rings(ring_list):
        clusters = []
        for ring in ring_list:
            merged = False
            for cluster in clusters:
                # if the current ring and cluster share at least 2 atoms consider them fused
                if len(cluster & ring) >= 2:
                    cluster |= ring
                    merged = True
                    break
            if not merged:
                clusters.append(set(ring))
        return clusters

    # Get all ring information from the molecule and store each ring as a set of atom indices.
    ring_info = mol.GetRingInfo()
    all_rings = [set(r) for r in ring_info.AtomRings()]
    if not all_rings:
        return False, "No ring systems detected; steroid nucleus expected in a sterol ester"

    # For each ester match, we try to see if its alcohol oxygen is attached (directly or via one extra bond)
    # to a candidate steroid nucleus.
    sterol_ester_found = False
    for match in ester_matches:
        # In our match, we expect the atoms to include one carbon (the carbonyl carbon) and one or more oxygens.
        # We now separate the candidate carbonyl carbon from the ester oxygens.
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
        # For each candidate oxygen, verify that its bond to the carbonyl carbon is a single bond.
        for oxy_idx in candidate_oxygens:
            bond = mol.GetBondBetweenAtoms(candidate_carbon, oxy_idx)
            if bond is None or bond.GetBondTypeAsDouble() != 1.0:
                continue  # skip if not a proper single bond
            # Now get the neighbors of this candidate oxygen (except the carbonyl carbon).
            o_atom = mol.GetAtomWithIdx(oxy_idx)
            non_carbonyl_neighbors = [nbr.GetIdx() for nbr in o_atom.GetNeighbors() if nbr.GetIdx() != candidate_carbon]
            if not non_carbonyl_neighbors:
                continue  # nothing attached to oxygen aside from the carbonyl
            # For each non-carbonyl neighbor, check if it is in at least one non-aromatic ring.
            for nbr_idx in non_carbonyl_neighbors:
                nbr_atom = mol.GetAtomWithIdx(nbr_idx)
                # Get rings that contain this neighbor atom.
                candidate_rings = []
                for ring in all_rings:
                    if nbr_idx in ring:
                        # Check that the ring is non-aromatic.
                        # We require that not every atom in the ring is aromatic.
                        # (Steroid rings are typically non-aromatic.)
                        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                        if not all(a.GetIsAromatic() for a in ring_atoms):
                            candidate_rings.append(ring)
                if not candidate_rings:
                    continue  # this neighbor is not in any non-aromatic ring so skip
                # Cluster the candidate rings for this neighbor.
                clusters = cluster_rings(candidate_rings)
                for cluster in clusters:
                    # Count how many rings from candidate_rings are fused into this cluster.
                    # We also check whether at least one ring in the cluster is 5-membered.
                    count = 0
                    has_five = False
                    for ring in candidate_rings:
                        if len(ring & cluster) >= 2:
                            count += 1
                            if len(ring) == 5:
                                has_five = True
                    # We now use a relaxed version of steroid nucleus criteria:
                    # at least 4 rings and one is five-membered.
                    if count >= 4 and has_five:
                        sterol_ester_found = True
                        break
                if sterol_ester_found:
                    break
            if sterol_ester_found:
                break
        if sterol_ester_found:
            break

    if not sterol_ester_found:
        return False, ("Ester group found but no candidate alcohol oxygen appears to be connected (directly or via one extra bond) "
                       "to a fused non-aromatic ring cluster having at least 4 rings (with at least one 5-membered ring) "
                       "as expected for a sterol nucleus")
    
    return True, ("Contains an ester group (C(=O)O) where the alcohol oxygen is connected directly or through a short bond "
                  "to a fused non-aromatic ring cluster with at least 4 rings (including at least one 5-membered ring), "
                  "consistent with a sterol ester")

# Example usage. You can test on examples such as:
if __name__ == '__main__':
    # Try an example SMILES (e.g., cholesteryl 11-hydroperoxy-eicosatetraenoate)
    smiles_example = "C(\\CCCC(=O)O[C@H]1CC[C@]2(C(C1)=CC[C@@]3([C@@]2(CC[C@]4([C@]3(CC[C@@]4([C@@](CCCC(C)C)(C)[H])[H])[H])C)[H])[H])C)=C\\C/C=C\\CCCCCOO"
    result, reason = is_sterol_ester(smiles_example)
    print("Result:", result)
    print("Reason:", reason)