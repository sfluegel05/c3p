"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Cyclohexenones
Definition: Any six‐membered alicyclic ketone having one C=C double bond (alkene) in the ring,
with a ketone functionality (C=O) placed on a ring carbon that is directly adjacent (in the ring)
to one of the alkene carbons (ensuring conjugation). In our new approach we only accept rings that 
are isolated (none of the atoms are shared with another six‐membered ring).
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES string.
      2. Get the ring information (using GetRingInfo).
      3. Loop over each six‐membered ring.
         a. Check that every atom in the ring is a carbon and non‐aromatic.
         b. Also require that each of the atoms in this ring is found in exactly one six‐membered ring
            (i.e. the ring is not fused with another six‐membered ring).
         c. Count the number of C=C double bonds that occur between two ring atoms (ignoring any C=O).
            Exactly one such alkene bond is required.
         d. For each ring atom, examine its neighbors that lie outside the ring. If an oxygen (atomic number 8)
            is double‐bonded (ketone) to a ring carbon then note that ring atom as “ketone‐bearing.”
         e. Finally, check “conjugation”: at least one ketone‐bearing ring atom must be adjacent (in the ring)
            to one of the carbons that forms the unique alkene bond.
      4. If any ring passes all checks, return True with an explanation.
    Returns:
      (bool, str): (True, explanation) if a valid cyclohexenone ring is detected; otherwise (False, explanation).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in molecule"

    # Helper: For each atom idx, count how many 6-membered rings (from our all_rings) it belongs to.
    ring_membership_count = {}
    for ring in all_rings:
        if len(ring) == 6:
            for idx in ring:
                ring_membership_count[idx] = ring_membership_count.get(idx, 0) + 1

    # Process each six-membered ring.
    for ring in all_rings:
        if len(ring) != 6:
            continue  # Only consider six-membered rings
            
        ring_set = set(ring)
        # Check that every atom is carbon, non-aromatic, and belongs to exactly one 6-membered ring.
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                valid_ring = False
                break
            if atom.GetIsAromatic():
                valid_ring = False
                break
            # For fused rings: if the atom belongs to >1 six-membered ring, skip.
            count = ring_membership_count.get(idx, 0)
            if count != 1:
                valid_ring = False
                break
        if not valid_ring:
            continue
            
        # Identify ring C=C bonds (alkene bonds). Only consider bonds where both atoms are in the ring.
        alkene_bonds = []
        for i in range(len(ring)):
            # To get proper ring connectivity, consider the ring as cyclic:
            a_idx = ring[i]
            b_idx = ring[(i+1) % len(ring)]  # next atom in the cycle
            bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is None:
                continue
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Make sure this is not a carbonyl bond: both atoms must be carbons (which they are in our check)
                alkene_bonds.append((a_idx, b_idx))
        # We require exactly one C=C double bond in the ring.
        if len(alkene_bonds) != 1:
            continue
        alkene_atoms = set(alkene_bonds[0])  # the two atoms forming the alkene

        # Identify ketone functionality. For each ring atom, check for an exocyclic oxygen that is double-bonded.
        ketone_atoms = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider substituents not in the ring.
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check that the bond is a double bond (C=O).
                    bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        ketone_atoms.add(idx)
        if not ketone_atoms:
            continue  # no ketone functionality found on this ring

        # Check conjugation:
        # We require that at least one ketone-carbon (in the ring) is directly bonded (in the ring) to one of the alkene carbons.
        conjugated = False
        for ket_idx in ketone_atoms:
            atom = mol.GetAtomWithIdx(ket_idx)
            # Find ring neighbors of ket_idx by intersecting its neighbors with ring_set.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set and nbr.GetIdx() in alkene_atoms:
                    conjugated = True
                    break
            if conjugated:
                break

        if conjugated:
            return True, ("Found an isolated six-membered non‐aromatic carbon ring having exactly one ring alkene "
                          "bond and a conjugated exocyclic ketone group; cyclohexenone detected.")
    return False, "No isolated cyclohexenone ring (six-membered ring with one ring alkene and a conjugated ketone) was found."

# Example test (uncomment to try):
# test_smiles = "O[C@H]1[C@H](O)C(=CC(=O)C1N)CO"  # 2-Aminovalienone example (true positive)
# result, reason = is_cyclohexenones(test_smiles)
# print(result, reason)