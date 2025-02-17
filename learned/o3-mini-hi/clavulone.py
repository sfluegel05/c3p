"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone – a class of esterified prostanoids obtained from marine corals.
Heuristic criteria (improved):
  1. The molecule must contain at least one 5‑membered ring (candidate cyclopentenone) that has:
       • exactly one exocyclic carbonyl group (an oxygen double‐bonded to a ring atom; outside the ring),
       • exactly one double bond connecting two ring atoms,
       • and the candidate ring should not be fused with a second ring.
       • Also, the ring atom bearing the exocyclic carbonyl is required to participate in the ring’s double bond (suggesting conjugation).
  2. The molecule must contain at least one ester substituent. Here, an ester is defined by the substructure
       “[#6]-[OX2]-[CX3](=O)[#6]”. Furthermore, at least one such ester must be “close” (graph distance ≤2 bonds)
       to the candidate ring.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_clavulone(smiles: str):
    """
    Determines if a molecule qualifies as a clavulone (esterified prostanoid) based on its SMILES string.

    The criteria (heuristic) are:
      1. The molecule contains at least one 5‐membered ring that:
           - has exactly one exocyclic carbonyl (an oxygen double‐bonded to a ring atom, where the O is not in the ring),
           - has exactly one double bond between atoms in the ring,
           - and the atom that bears the exocyclic carbonyl also is part of that internal double bond.
           - In addition, this candidate ring should not be fused with another ring.
      2. The molecule contains at least one ester substituent – defined as a substructure matching
         “[#6]-[OX2]-[CX3](=O)[#6]” – that is “close” (graph distance ≤2) to one of the candidate ring atoms.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule satisfies our clavulone heuristic.
        str : Explanation for the decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_ring = None
    candidate_ring_idx = None
    # Loop over 5-membered rings and check for a cyclopentenone-type candidate.
    for ring in all_rings:
        if len(ring) != 5:
            continue

        # 1a. Count exocyclic carbonyl groups for atoms in ring.
        exo_carbonyl_count = 0
        # Record which ring atoms bear an exocyclic carbonyl.
        carbonyl_bearing_atoms = []
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check for oxygen attached by a double bond.
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        exo_carbonyl_count += 1
                        carbonyl_bearing_atoms.append(atom_idx)
        if exo_carbonyl_count != 1:
            continue

        # 1b. Count internal (ring–ring) double bonds.
        ring_atom_set = set(ring)
        internal_double_count = 0
        # Also record bonds (as pairs of indices) that are double.
        double_bonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    internal_double_count += 1
                    double_bonds.append((a1, a2))
        if internal_double_count != 1:
            continue

        # 1c. Check for conjugation:
        # The ring atom bearing the exocyclic carbonyl must be one end of the internal double bond.
        conjugated = False
        for db in double_bonds:
            if carbonyl_bearing_atoms[0] in db:
                conjugated = True
                break
        if not conjugated:
            continue

        # 1d. Check that the candidate ring is not fused with any other ring.
        fused = False
        for other_ring in all_rings:
            if other_ring == ring:
                continue
            # If the intersection between rings has more than one atom, consider it fused.
            if len(set(ring).intersection(set(other_ring))) > 1:
                fused = True
                break
        if fused:
            continue

        # This ring meets our cyclopentenone-type criteria.
        candidate_ring = ring
        candidate_ring_idx = set(ring)
        break

    if candidate_ring is None:
        msg = ("No suitable cyclopentenone-type ring found: need a 5-membered ring with exactly one exocyclic carbonyl "
               "and one internal double bond (with conjugation) and not fused with another ring.")
        return False, msg

    # 2. Look for at least one ester substituent.
    # Define an ester SMARTS pattern:
    ester_smarts = "[#6]-[OX2]-[CX3](=O)[#6]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    if ester_query is None:
        return False, "Error generating ester SMARTS pattern"
    
    ester_matches = mol.GetSubstructMatches(ester_query)
    ester_found = False
    # For each ester match, check if any atom in the match is within graph distance 2 of any candidate ring atom.
    for match in ester_matches:
        # match is a tuple of atom indices corresponding to the pattern.
        # We want to check, for example, the oxygen in the ester (pattern position 1)
        for atom_match in match:
            for ring_atom in candidate_ring:
                # Get the shortest path (list of atom indices) between ring_atom and the ester atom.
                path = rdmolops.GetShortestPath(mol, ring_atom, atom_match)
                # The number of bonds equals (len(path)-1)
                if len(path) - 1 <= 2:
                    ester_found = True
                    break
            if ester_found:
                break
        if ester_found:
            break

    if not ester_found:
        msg = ("Candidate cyclopentenone-type ring found but no ester substituent found in close proximity "
               "(graph distance ≤2) to the ring.")
        return False, msg

    # If we arrived here, the molecule meets our heuristic for clavulone.
    msg = ("Contains a cyclopentenone-type 5-membered ring (isolated and conjugated to an exocyclic carbonyl) "
           "plus an ester substituent in close proximity – consistent with clavulone.")
    return True, msg


# Example usage (for testing):
if __name__ == "__main__":
    # Testing with one of the clavulone examples: punaglandin 2
    smiles_example = ("ClC=1C(=O)[C@@]([C@@](O)(C/C=C\\CCCCC)C1)"
                      "([C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O)[H]")
    result, reason = is_clavulone(smiles_example)
    print(result, reason)