"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugates

A beta-D-glucosiduronic acid (glucuronide) conjugate is defined as a glucuronic acid unit 
that has undergone formal condensation with another substance via its anomeric position.
Rather than using a single rigid SMARTS pattern, this implementation first searches for 
a pyranose-like ring (six-membered ring with one oxygen) that possesses a carboxyl substituent 
(evidence of oxidation to glucuronic acid) and then checks that one of its ring carbons has an 
exocyclic oxygen (the anomeric oxygen) that is linked to an external (non-sugar) heavy fragment.
Note that this heuristic approach may not capture every nuance and might still miss edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid conjugate (glucuronide) based on its SMILES string.
    
    Approach:
      1. Parse the molecule and add explicit hydrogens.
      2. Identify candidate pyranose rings: rings of 6 atoms that include exactly one oxygen.
      3. For each candidate ring:
         a. Look for a ring carbon that bears an exocyclic carboxyl group. This group should be seen as a
            carbon not in the ring that has one double-bonded O and one single-bonded O (which might be deprotonated).
         b. Look for an "anomeric" oxygen: an exocyclic oxygen (attached to a ring carbon) that is not part of the ring.
            Check that the substituent linked to that oxygen (other than the ring carbon) has heavy‐atomic (non‐H)
            connectivity. (We also reject common phosphate linkers.)
      4. If a candidate ring meets both conditions, classify the molecule as a beta-D-glucosiduronic acid conjugate.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: A textual reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to clarify attachment.
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    found_candidate = False
    candidate_reason = ""
    
    # helper: check if atom is a carboxyl carbon (C(=O)[O,OH])
    def is_carboxyl(carbon):
        # Must be carbon
        if carbon.GetAtomicNum() != 6:
            return False
        double_O = 0
        single_O = 0
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    double_O += 1
                elif bond.GetBondType() == rdchem.BondType.SINGLE:
                    single_O += 1
        return (double_O >= 1 and single_O >= 1)
    
    # For each ring, check if it is a candidate sugar ring.
    for ring in rings:
        if len(ring) != 6:
            continue  # we require a 6-membered ring
        # Count the oxygens in the ring.
        ring_oxygens = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue  # typical pyranose should have one ring oxygen
        # Now identify a carboxyl substituent (exocyclic group) attached to one of the ring carbons.
        glucuronic_found = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Consider only carbons from the ring.
            if atom.GetAtomicNum() != 6:
                continue
            # Look at neighbors not in ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # A carboxyl group is usually represented on a carbon substituent.
                if nbr.GetAtomicNum() == 6 and is_carboxyl(nbr):
                    glucuronic_found = True
                    break
            if glucuronic_found:
                break
        if not glucuronic_found:
            # This ring does not look like glucuronic acid (missing carboxyl group).
            continue

        # Next, check for evidence of conjugation via the anomeric oxygen.
        # Search among ring carbons for one that has an exocyclic oxygen.
        anomeric_found = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look among its neighbors for an oxygen not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() != 8:
                    continue
                # Likely candidate for anomeric oxygen.
                # Now, check if this oxygen is attached to something outside the sugar ring.
                external_neighbors = []
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == atom.GetIdx():
                        continue
                    # Ignore hydrogens.
                    if nbr2.GetAtomicNum() == 1:
                        continue
                    external_neighbors.append(nbr2)
                # If the oxygen links to an external group with at least one heavy neighbor, we take that as evidence of conjugation.
                # Also, if the external neighbor is a phosphate (P) atom, then skip because that may be part of a UDP-sugar.
                for ext in external_neighbors:
                    if ext.GetAtomicNum() == 15:  # phosphorus
                        continue
                    # Check that ext itself is not just a lone oxygen (should be a heavier fragment)
                    heavy_count = sum(1 for n in ext.GetNeighbors() if n.GetAtomicNum() > 1)
                    if heavy_count >= 1:
                        anomeric_found = True
                        break
                if anomeric_found:
                    break
            if anomeric_found:
                break

        if anomeric_found:
            found_candidate = True
            candidate_reason = ("Found candidate sugar ring with glucuronic acid (carboxyl group detected) "
                                "and anomeric oxygen conjugated to an external substituent with sufficient heavy connectivity.")
            break  # if one candidate ring qualifies, we classify as a glucuronide conjugate

    if found_candidate:
        return True, candidate_reason
    else:
        # We distinguish here whether we found any glucuronic‐sugar ring at all.
        for ring in rings:
            if len(ring) == 6:
                ring_oxygens = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]
                if len(ring_oxygens) == 1:
                    # If such a ring was found but no conjugation evidence was present…
                    return False, ("Beta-D-glucuronic acid-like ring detected but no external substituent from the anomeric oxygen "
                                   "was found. Likely the sugar is present free rather than conjugated.")
        return False, "No beta-D-glucuronic acid moiety found in the molecule"

    
# Example usage when run as a script.
if __name__ == "__main__":
    # Test with one of the provided examples: tamoxifen N-beta-D-glucosiduronic acid
    test_smiles = "CC\\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1"
    result, reason = is_beta_D_glucosiduronic_acid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)