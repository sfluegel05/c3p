"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: Polysaccharide
Definition: A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
            For our purposes, we assume a polysaccharide has (roughly) more than 10 sugar rings.
The function is_polysaccharide takes a SMILES string as input and returns a tuple (bool, reason).
This version improves by relaxing strict atom‐criteria and by checking that candidate rings have
a reasonable number of hydroxyl substituents.
"""

from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.

    For our improved heuristic, we first extract all rings. Then among rings that are of size 5 or 6,
    we count as “sugar‐like” rings if:
      • The ring has at least one oxygen (and not too many extra O atoms, say no more than 2 in a 6‐membered ring)
      • It has at least 70% of its atoms as carbons (allowing for modifications)
      • At least two atoms in the ring bear a hydroxyl substituent (an –OH not in the ring)
    The molecule is classified as a polysaccharide if more than 10 (i.e. >= 11) such rings are found.
    (Note: In our test dataset, some true polysaccharides are borderline; one might adjust the threshold.)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a polysaccharide, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in molecule"

    sugar_ring_count = 0

    for ring in rings:
        ring_size = len(ring)
        # Only consider rings of size 5 or 6 (furanoses and pyranoses)
        if ring_size not in (5, 6):
            continue

        # Count ring oxygens and carbons and also check hydroxyl substituents outside the ring.
        ring_oxygen_count = 0
        ring_carbon_count = 0
        hydroxyl_substituents = 0

        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                ring_oxygen_count += 1
            elif atomic_num == 6:
                ring_carbon_count += 1
            # Optionally ignore other atoms

            # Check neighbors not in ring: if a neighbor is oxygen and is attached as a single bond (likely –OH)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring and neighbor.GetAtomicNum() == 8:
                    # Check that the bond is a single bond and that neighbor does not carry charges
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() == 1.0:
                        # To avoid counting carboxylates (the C=O oxygen), we require that the neighbor
                        # only has a single heavy-atom neighbor (the ring atom)
                        if len([a for a in neighbor.GetNeighbors() if a.GetAtomicNum() > 1]) == 1:
                            hydroxyl_substituents += 1

        # For a furanose ring (5-membered) one expects 1 oxygen and 4 carbons,
        # but we relax for slight modifications.
        if ring_size == 5:
            if ring_oxygen_count >= 1 and ring_carbon_count >= 3 and hydroxyl_substituents >= 1:
                sugar_ring_count += 1
        elif ring_size == 6:
            # For pyranoses, typically 1 oxygen and 5 carbons are expected.
            # Allow up to 2 oxygens in the ring (to tolerate minor substitutions) and require at least 4 carbons.
            if ring_oxygen_count >= 1 and ring_oxygen_count <= 2 and ring_carbon_count >= 4 and hydroxyl_substituents >= 2:
                sugar_ring_count += 1

    # We adjust the threshold a bit in case strict count misses borderline cases.
    threshold = 11  # more than 10 rings needed
    if sugar_ring_count >= threshold:
        return True, f"Contains {sugar_ring_count} sugar-like rings (monosaccharide residues) consistent with a polysaccharide."
    else:
        return False, f"Only {sugar_ring_count} sugar-like rings detected; need at least {threshold} for polysaccharide classification."


# Example usage:
if __name__ == "__main__":
    # Test with a typical glucose ring (1 unit)
    test_smiles = "OC1OC(O)C(O)C(O)C1O"  # one pyranose unit
    result, reason = is_polysaccharide(test_smiles)
    print("Is polysaccharide:", result)
    print("Reason:", reason)