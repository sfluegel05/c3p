"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Cyclohexenones
Definition: Any six‐membered alicyclic ketone having one C=C double bond (alkene) in the ring,
with a ketone functionality (C=O) placed on a ring carbon that is directly adjacent (in the ring)
to one of the alkene carbons (ensuring conjugation). We now allow cases where the cyclohexenone
substructure appears even if an atom is shared, provided that no bond (i.e. two adjacent ring atoms)
is shared between two six‐membered rings.
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule contains a cyclohexenone substructure based on its SMILES string.

    The algorithm:
      1. Parse the SMILES string.
      2. Get ring information (using GetRingInfo).
      3. Loop over each six‐membered ring:
         a. Ensure all atoms in the ring are carbons and non‐aromatic.
         b. Check if any bond in the ring is shared with another six‐membered ring (fused bond);
            if so, skip this ring.
         c. Count the number of C=C double bonds (alkene bonds) among ring–internal bonds.
            Exactly one such alkene bond is required.
         d. For each ring atom, check its neighbors outside the ring for an oxygen that is double‐bonded.
            Mark these ring atoms as “ketone‐bearing.”
         e. Check conjugation: at least one ketone–bearing ring atom must be adjacent (within the ring)
            to one of the two alkene atoms.
      4. If any ring passes, return True with an explanation; otherwise, return False.

    Returns:
      (bool, str): (True, explanation) if a valid cyclohexenone structure is detected; otherwise (False, explanation).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in molecule"

    # Gather only six-membered rings (list of tuples of atom indices)
    six_membered_rings = [ring for ring in all_rings if len(ring) == 6]
    if not six_membered_rings:
        return False, "No six-membered rings found in molecule"

    # For helping determine bond fusion, construct a set for each six-membered ring's bonds.
    # Each bond is represented as a frozenset({a_idx, b_idx}) for adjacent atoms in the ring.
    rings_with_bonds = []
    for ring in six_membered_rings:
        # sort the ring into a cyclic order using the ordering from GetRingInfo 
        rset = set(ring)
        bonds = set()
        n = len(ring)
        for i in range(n):
            a = ring[i]
            b = ring[(i+1)%n]
            bonds.add(frozenset((a, b)))
        rings_with_bonds.append((set(ring), bonds))

    # Now process each candidate six-membered ring.
    for ring in six_membered_rings:
        ring_set = set(ring)
        # (a) Check that every atom in the ring is carbon and non‐aromatic.
        atoms_ok = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                atoms_ok = False
                break
        if not atoms_ok:
            continue

        # (b) Check for bond fusion: for every adjacent pair in the ring, if that bond appears in any other six‐membered ring, consider the bond fused.
        # Flag if any bond is found in more than one ring.
        fused = False
        # First get bonds for this ring.
        bonds_this_ring = set()
        n = len(ring)
        for i in range(n):
            a = ring[i]
            b = ring[(i+1)%n]
            bonds_this_ring.add(frozenset((a, b)))
        # Now check against all other six-membered rings.
        for other_ring_set, other_bonds in rings_with_bonds:
            # Skip comparing the ring with itself
            if other_ring_set == ring_set:
                continue
            # If there is an edge (bond) common to both rings, consider it fused.
            if bonds_this_ring.intersection(other_bonds):
                fused = True
                break
        if fused:
            # Skip fused six-membered rings
            continue

        # (c) Identify ring alkene bonds.
        alkene_bonds = []
        for i in range(n):
            a_idx = ring[i]
            b_idx = ring[(i+1)%n]
            bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is None:
                continue
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Even though both atoms are carbons per our earlier check, we assume a C=O would not occur inside the ring.
                alkene_bonds.append((a_idx, b_idx))
        if len(alkene_bonds) != 1:
            continue

        # Save the two atom indices involved in the alkene bond
        alkene_atoms = set(alkene_bonds[0])

        # (d) Identify exocyclic ketone functionality.
        # For each ring atom, check for a neighbor not in the ring that is oxygen and double-bonded.
        ketone_atoms = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue  # only exocyclic groups
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        ketone_atoms.add(idx)
        if not ketone_atoms:
            continue  # no ketone functionality found on the ring

        # (e) Check for conjugation: at least one ketone-bearing ring atom must be adjacent (within the ring)
        # to one of the alkene atoms.
        conjugated = False
        for idx in ketone_atoms:
            # get ring neighbors for this index:
            atom_neighbors_in_ring = []
            for i in range(n):
                if ring[i] == idx:
                    # adjacent vertices in the cyclic order:
                    prev_atom = ring[(i-1)%n]
                    next_atom = ring[(i+1)%n]
                    atom_neighbors_in_ring.extend([prev_atom, next_atom])
            # Check if any neighbor is in the alkene_atoms
            for n_idx in atom_neighbors_in_ring:
                if n_idx in alkene_atoms:
                    conjugated = True
                    break
            if conjugated:
                break

        if conjugated:
            return True, ("Found a six‐membered, non‐aromatic carbon ring with exactly one ring alkene bond and an exocyclic ketone "
                          "group conjugated to the alkene; cyclohexenone substructure detected.")
    return False, "No valid cyclohexenone ring (six‐membered ring with one alkene bond and a conjugated ketone) was found."

# Example test (uncomment to try):
# test_smiles = "O[C@H]1[C@H](O)C(=CC(=O)C1N)CO"  # 2-Aminovalienone example: expected True
# result, reason = is_cyclohexenones(test_smiles)
# print(result, reason)