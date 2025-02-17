"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: Sterol Ester
Definition: A steroid ester obtained by formal condensation of the carboxy group 
of any carboxylic acid with the 3-hydroxy group of a sterol.

Heuristic:
  1. The molecule must contain an ester functional group (C(=O)O).
  2. The molecule must have a steroid nucleus – here at least one 5-membered ring 
     and three 6-membered rings.
  3. For each ester group, identify the “ester oxygen” (the oxygen attached through a single bond
     to the carbonyl carbon). Then, require that one of its neighbors (other than the carbonyl carbon)
     is in one of the rings that is part of a fused (i.e. multiple-ring) system.
"""

from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    
    A sterol ester is defined as a steroid ester formed by reaction of any 
    carboxylic acid with the 3-hydroxy group of a sterol. We expect an ester (C(=O)O)
    present and that the ester oxygen (from the alcohol side) is attached directly or 
    via a short bond to a fused ring system characteristic of a steroid nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sterol ester, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of an ester group using a SMARTS pattern: C(=O)O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group (C(=O)O) found"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Count rings by size to (roughly) ensure a steroid nucleus.
    num_5_membered = sum(1 for ring in rings if len(ring) == 5)
    num_6_membered = sum(1 for ring in rings if len(ring) == 6)
    if num_5_membered < 1 or num_6_membered < 3:
        return False, ("Steroid nucleus not detected – expected at least one 5-membered ring and three 6-membered rings, "
                       f"found {num_5_membered} and {num_6_membered}")
    
    # Build a mapping for each atom: count how many 5- or 6-membered rings the atom is in.
    ring_counts = {i: 0 for i in range(mol.GetNumAtoms())}
    for ring in rings:
        if len(ring) in (5,6):
            for idx in ring:
                ring_counts[idx] += 1
    
    # Now, among the ester groups we search for one that comes from esterification
    # at a sterol (3-hydroxy) position.
    # For each ester match, do:
    #   - Identify the carbonyl carbon (by atomic number 6) and the oxygen atoms.
    #   - Choose an oxygen atom that is attached by a single bond to the carbonyl carbon.
    #   - Then, check if one of its neighbors (aside from the carbonyl carbon) belongs
    #     to a fused ring (i.e. ring_counts at least 2).
    sterol_ester_found = False
    for match in ester_matches:
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
        
        # For each oxygen attached to the carbonyl carbon, look for a single bond.
        for oxy_idx in candidate_oxygens:
            bond = mol.GetBondBetweenAtoms(candidate_c, oxy_idx)
            if bond is None:
                continue
            # Accept only the oxygen that is single-bonded to the carbonyl carbon.
            if bond.GetBondTypeAsDouble() != 1.0:
                continue
            ester_oxygen_idx = oxy_idx
            ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)
            # Get neighbors of the ester oxygen excluding the carbonyl carbon.
            neighbor_idxs = [nbr.GetIdx() for nbr in ester_oxygen.GetNeighbors() if nbr.GetIdx() != candidate_c]
            if not neighbor_idxs:
                continue
            # Now check if any neighbor is in a fused ring by requiring that
            # the neighbor is in a 5- or 6-membered ring and appears in at least 2 such rings.
            for nbr_idx in neighbor_idxs:
                if ring_counts.get(nbr_idx, 0) >= 2:
                    sterol_ester_found = True
                    break
            if sterol_ester_found:
                break
        if sterol_ester_found:
            break
    
    if not sterol_ester_found:
        return False, ("Ester group found but the candidate ester oxygen does not appear "
                       "to connect to a fused steroid ring system (expected for a sterol-derived alcohol)")
    
    return True, ("Contains an ester group with an alcohol side linked (directly or via a short bond) "
                  "to a fused steroid nucleus, consistent with a sterol ester")
    
# Example usage:
if __name__ == "__main__":
    # Example: use one of the provided test SMILES (for instance, CE(20:2(6Z,9Z)))
    smiles_example = "O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O"
    result, reason = is_sterol_ester(smiles_example)
    print("Result:", result)
    print("Reason:", reason)