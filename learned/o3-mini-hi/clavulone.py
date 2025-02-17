"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone – a class of esterified prostanoids obtained from marine corals.
Heuristic criteria:
  1. The molecule must contain a 5–membered (cyclopentenone‐type) ring that has:
       - exactly one exocyclic carbonyl group (an oxygen double‐bonded to a ring atom, but not within the ring)
       - exactly one double bond connecting ring atoms.
  2. At least one ester substituent – defined as an oxygen (bonded by a single bond to a ring atom) 
     that is further attached to a carbon which itself is double‐bonded to an oxygen – must be attached 
     to one of the ring atoms.
"""

from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule qualifies as a clavulone (esterified prostanoid) based on its SMILES string.

    The criteria are:
      1. The molecule contains at least one 5‐membered ring such that:
           - exactly one exocyclic carbonyl is attached to a ring atom (an oxygen double‐bonded 
             to a ring atom from outside the ring)
           - exactly one double bond exists between ring atoms.
      2. At least one ester substituent is found attached to one of the ring atoms. Here an ester substituent
         is defined as an oxygen (attached by a single bond to the ring) that is connected to a carbon which,
         in turn, is double‐bonded to an oxygen.

    Args:
        smiles (str): SMILES string for the molecule

    Returns:
        bool: True if the molecule satisfies our clavulone criteria, False otherwise.
        str : Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_ring = None

    # Loop over all rings to find a candidate which is 5-membered
    for ring in all_rings:
        if len(ring) != 5:
            continue

        # 1a. Count exocyclic carbonyls for atoms in ring.
        exo_carbonyl_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # look at neighbors not in ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # if the neighbor is oxygen and attached by a double bond (C=O)
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        exo_carbonyl_count += 1
        
        if exo_carbonyl_count != 1:
            continue

        # 1b. Count internal (ring–ring) double bonds.
        ring_atom_set = set(ring)
        internal_double_count = 0
        # Only consider bonds where both atoms are in the ring:
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    internal_double_count += 1
        if internal_double_count != 1:
            continue

        # Found a candidate ring that seems cyclopentenone‐like.
        candidate_ring = ring
        break

    if candidate_ring is None:
        msg = ("No suitable cyclopentenone-type ring found: "
               "need a 5-membered ring with exactly one exocyclic carbonyl "
               "and one internal double bond.")
        return False, msg

    # 2. Look for at least one ester substituent attached to any candidate ring atom.
    ester_found = False
    for atom_idx in candidate_ring:
        ring_atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in ring_atom.GetNeighbors():
            # Only consider neighbors not part of the candidate ring.
            if nbr.GetIdx() in candidate_ring:
                continue
            # The ester oxygen is expected to be an oxygen atom, connected by a single bond.
            if nbr.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
            if not bond or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Check the connectivity of this oxygen:
            # For an ester substituent, the oxygen should have exactly 2 neighbors (one is the ring atom).
            if nbr.GetDegree() < 2:
                continue
            # Find the neighbor (other than our ring atom) attached to this oxygen.
            for oxy_nbr in nbr.GetNeighbors():
                if oxy_nbr.GetIdx() == atom_idx:
                    continue
                # That neighbor should be a carbon (the carbonyl carbon of the ester).
                if oxy_nbr.GetAtomicNum() != 6:
                    continue
                bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), oxy_nbr.GetIdx())
                if not bond2 or bond2.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Check that this carbon is double-bonded to an oxygen (i.e. forming a carbonyl).
                has_carbonyl = False
                for c_nbr in oxy_nbr.GetNeighbors():
                    if c_nbr.GetIdx() == nbr.GetIdx():
                        continue
                    if c_nbr.GetAtomicNum() != 8:
                        continue
                    bond3 = mol.GetBondBetweenAtoms(oxy_nbr.GetIdx(), c_nbr.GetIdx())
                    if bond3 and bond3.GetBondType() == Chem.BondType.DOUBLE:
                        has_carbonyl = True
                        break
                if has_carbonyl:
                    ester_found = True
                    break
            if ester_found:
                break
        if ester_found:
            break

    if not ester_found:
        msg = "Candidate cyclopentenone-type ring found but no ester substituent attached to it."
        return False, msg

    # If we reached here, we have met our heuristic.
    return True, ("Contains a cyclopentenone-type 5-membered ring with one exocyclic carbonyl "
                   "and one internal double bond plus an ester substituent attached to the ring – "
                   "consistent with clavulone.")

# Example usage (for testing):
if __name__ == "__main__":
    # Here is one test with punaglandin 2 (a known clavulone example).
    smiles_example = "ClC=1C(=O)[C@@]([C@@](O)(C/C=C\\CCCCC)C1)([C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O)[H]"
    result, reason = is_clavulone(smiles_example)
    print(result, reason)