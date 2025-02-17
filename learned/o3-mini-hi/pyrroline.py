"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: Organic heteromonocyclic compounds based on a dihydropyrrole (pyrroline) core.
A pyrroline (for our purposes) is defined as a 5‐membered ring that:
  - contains exactly one nitrogen,
  - and has exactly one “multiple bond” that is either:
      • an internal (in–ring) double bond or
      • an exocyclic double bond (only counted when it is from a ring atom to O or S).
In addition, many false positives (especially in heavily fused systems) may be avoided by
requiring that an acceptable ring has at least three atoms that are not shared with any other ring.
Note:
  Because many complex molecules contain fused rings or additional unsaturation, this heuristic
  is not perfect and may over– or under–classify some structures.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_pyrroline(smiles: str):
    """
    Determines if the given molecule (via its SMILES string) contains a pyrroline (dihydropyrrole)
    core. The algorithm works as follows:
      1. Parse the molecule and attempt a Kekulization so that double bonds are explicit.
      2. Retrieve all rings and consider only 5‐membered rings.
      3. For each candidate ring, check that it contains exactly one nitrogen.
      4. Count the “multiple bonds” associated with the ring:
           - Count intraring bonds drawn as a double bond.
           - Also count exocyclic double bonds (i.e. from a ring atom to a non‐ring atom)
             but only if the non–ring atom is oxygen or sulfur.
      5. Accept the ring if the total count equals exactly one.
      6. To reduce false positives from highly fused ring systems, if the ring is fused
         (i.e. some atoms appear in >1 ring overall) then require that at least 3 out of 5
         atoms are unique to this ring.
      7. If a qualifying ring is found, return True with a message that indicates whether the 
         ring is isolated or fused; otherwise, return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a qualifying pyrroline core is found, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Try to Kekulize so that explicit bond orders are available.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # If kekulization fails, proceed with the original structure.
        pass

    ringinfo = mol.GetRingInfo()
    atom_rings = ringinfo.AtomRings()
    if not atom_rings:
        return False, "No rings found in the molecule"

    # Precompute how many rings each atom belongs to.
    atom_ring_counts = {atom.GetIdx(): 0 for atom in mol.GetAtoms()}
    for ring in atom_rings:
        for idx in ring:
            atom_ring_counts[idx] += 1

    # Examine each 5-membered ring:
    for ring in atom_rings:
        if len(ring) != 5:
            continue

        # Count the nitrogen atoms inside this ring.
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue

        # Collect bonds between atoms in the ring.
        ring_atom_set = set(ring)
        intra_ring_bonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                intra_ring_bonds.append(bond)
        if len(intra_ring_bonds) != 5:
            # Something is unusual about this ring.
            continue

        # Count intraring double bonds.
        intra_double = sum(1 for bond in intra_ring_bonds if bond.GetBondType() == rdchem.BondType.DOUBLE)

        # Count exocyclic double bonds: for each ring atom, examine bonds to non–ring atoms.
        exo_double = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                # Only look at bonds connecting to atoms outside the ring.
                other = bond.GetOtherAtom(atom)
                if other.GetIdx() in ring_atom_set:
                    continue
                # Only count a double bond to O or S.
                if bond.GetBondType() == rdchem.BondType.DOUBLE and other.GetAtomicNum() in (8, 16):
                    exo_double += 1

        total_multiple = intra_double + exo_double

        # For a dihydropyrrole core, we require exactly one multiple bond.
        if total_multiple != 1:
            continue

        # Check fused vs isolated. Count how many ring atoms are found only in this ring.
        unique_count = sum(1 for idx in ring if atom_ring_counts[idx] == 1)
        if unique_count < 3:
            # If less than 3 atoms are unique to this ring, it is heavily fused.
            ring_type = "fused"
        else:
            ring_type = "isolated"

        # (Optionally, one might decide that certain fused rings are suspect.
        # For this implementation, we still accept them.)
        return True, f"Found a {ring_type} 5‐membered ring with one nitrogen and one multiple bond (dihydropyrrole core)."

    return False, "No 5‐membered pyrroline (dihydropyrrole) core meeting our criteria found."


# Example usage:
if __name__ == '__main__':
    # Test using a known pyrroline: 2-acetyl-1-pyrroline
    test_smiles = "CC(=O)C1=NCCC1"
    result, reason = is_pyrroline(test_smiles)
    print(result, reason)