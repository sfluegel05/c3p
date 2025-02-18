"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: Any six-membered alicyclic ketone having one double bond in the ring (cyclohexenones).

A cyclohexenone motif in our definition is a non-aromatic six-membered carbon ring that
contains exactly one C=C double bond and one of the ring carbons bears an exocyclic 
carbonyl group (C=O) that is directly conjugated to that double bond.
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.

    A cyclohexenone is defined as any six-membered alicyclic ketone having one double
    bond in the ring and a C=O group conjugated (directly attached) to the alkene.

    The algorithm will:
      - Parse the SMILES.
      - Loop over all rings found.
      - For rings of exactly six atoms:
          * Ensure every atom in the ring is a non-aromatic carbon.
          * Identify all bonds within the ring that are double bonds between carbons.
            We require exactly one such ring double bond.
          * Look for a carbonyl group at a ring carbon (a double bond from a ring carbon
            to an oxygen that is not part of the ring).
          * Require that the carbonyl-bearing carbon is one of the two atoms of the ring alkene.
          * If such a ring is found, classify as cyclohexenone.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        bool: True if the molecule is classified as cyclohexenone, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule"

    # Loop over each ring
    for ring in ring_info:
        if len(ring) != 6:
            continue  # only consider six-membered rings

        # Ensure all ring atoms are carbon and non-aromatic
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

        # Determine the bonds that lie completely inside the ring.
        ring_set = set(ring)
        ring_double_bonds = []  # list of tuples: (bond, idx1, idx2)
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                # Count only double bonds between carbons
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    atom1 = mol.GetAtomWithIdx(a1)
                    atom2 = mol.GetAtomWithIdx(a2)
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        ring_double_bonds.append((bond, a1, a2))
        # We expect exactly one carbon-carbon double bond in the ring.
        if len(ring_double_bonds) != 1:
            continue

        # For each ring atom, check if it has a carbonyl group (C=O outside the ring).
        # Record the indices of ring atoms that have a carbonyl double bond.
        carbonyl_atoms = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Skip if the neighbor is in the ring.
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Get the bond between the ring atom and the oxygen.
                    bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        carbonyl_atoms.add(idx)
        if not carbonyl_atoms:
            continue  # no carbonyl found on any ring atom

        # Check conjugation: One of the carbonyl atoms must be part of the ring alkene.
        _, alkene_idx1, alkene_idx2 = ring_double_bonds[0]
        if alkene_idx1 in carbonyl_atoms or alkene_idx2 in carbonyl_atoms:
            return True, ("Found a six-membered non-aromatic carbon ring with one C=C bond and a "
                          "conjugated ketone (C=O) attached to an alkene carbon; cyclohexenone motif detected.")
        # If no carbonyl-bearing atom is part of the alkene, try next ring.
    return False, "No six-membered ring with a conjugated alkene and ketone (cyclohexenone) was found."

# Example test calls (uncomment to run)
# test_smiles = "O[C@H]1[C@H](O)C(=CC(=O)C1N)CO"  # 2-Aminovalienone example
# result, reason = is_cyclohexenones(test_smiles)
# print(result, reason)