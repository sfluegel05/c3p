"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Any six-membered alicyclic ketone having one double bond in the ring (cyclohexenones).

A cyclohexenone motif in our definition is a non-aromatic six-membered carbon ring that:
  - Contains exactly one carbon-carbon double bond (alkene) inside the ring.
  - Contains a ketone group where a ring carbon is double-bonded to an oxygen (not part of the ring).
  - The carbonyl-bearing ring carbon must be directly adjacent (in the ring) to one of the carbons 
    involved in the alkene double bond (ensuring conjugation).
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.

    A cyclohexenone is defined as any six-membered alicyclic ketone having one C=C double bond in 
    the ring and a carbonyl group (C=O) placed on a ring carbon conjugated with the alkene.

    The algorithm:
      1. Parses the SMILES.
      2. Iterates over all rings (using GetRingInfo()).
      3. For each six-membered ring composed entirely of non-aromatic carbons:
         a. Count the number of C=C double bonds inside the ring (ignore the C=O bond).
         b. Identify ring atoms that have a double bond to oxygen with the O not being in the ring 
            (the ketone functionality).
         c. Verify that one of these ketone-bearing atoms is adjacent (within the ring) to one of the 
            carbons forming the alkene double bond.
    Returns:
      bool: True if classified as cyclohexenone, False otherwise.
      str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from molecule.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in molecule"

    # Loop over each ring
    for ring in rings:
        if len(ring) != 6:
            continue  # Only consider six-membered rings

        # Ensure that every atom in the ring is a carbon and the ring is non-aromatic.
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                valid_ring = False
                break
            if atom.GetIsAromatic():
                valid_ring = False
                break
        if not valid_ring:
            continue

        ring_set = set(ring)
        # Identify ring C=C bonds (alkene bonds, within the ring).
        alkene_bonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Both atoms must be carbon for a C=C alkene.
                    if (mol.GetAtomWithIdx(a1).GetAtomicNum() == 6 and
                        mol.GetAtomWithIdx(a2).GetAtomicNum() == 6):
                        alkene_bonds.append((a1, a2))
        # We require exactly one C=C double bond within the ring.
        if len(alkene_bonds) != 1:
            continue

        alkene_atoms = set(alkene_bonds[0])  # the two atoms forming the alkene

        # Identify a ketone functionality on the ring.
        # For each ring atom, check neighbors not in the ring: if an oxygen is double-bonded.
        ketone_atoms = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue  # only consider substituents outside the ring
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        ketone_atoms.add(idx)
        if not ketone_atoms:
            continue  # This ring does not have a ketone functional group

        # Check conjugation: the ketone-bearing ring carbon must be adjacent (inside the ring)
        # to one of the alkene carbons.
        conjugated = False
        for ket_idx in ketone_atoms:
            atom = mol.GetAtomWithIdx(ket_idx)
            # Get the atoms in the ring that are connected to ket_idx:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set and nbr.GetIdx() in alkene_atoms:
                    conjugated = True
                    break
            if conjugated:
                break

        if conjugated:
            return True, "Found six-membered non-aromatic carbon ring with one ring alkene and a conjugated ketone group; cyclohexenone detected."
    return False, "No six-membered ring with a conjugated alkene and ketone (cyclohexenone) was found."

# Example test (uncomment to try):
# test_smiles = "O[C@H]1[C@H](O)C(=CC(=O)C1N)CO"  # 2-Aminovalienone example
# result, reason = is_cyclohexenones(test_smiles)
# print(result, reason)