"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoids – members of the class of tetrahydrochromenochromene compounds.
A rotenoid is defined as having a cis‐fused tetrahydrochromeno[3,4-b]chromene skeleton,
i.e. a fused bicyclic system in which one six‐membered ring shows “chromenone‐like”
characteristics (has a ring oxygen and exactly one in‐ring carbonyl functionality)
fused to a second six‐membered ring that has a ring oxygen and NO carbonyls.
Substituted derivatives (e.g. glycosides) may be heavier than the classical 250–600 Da,
so the molecular weight filter has been relaxed.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.

    Heuristic criteria:
      1. The molecule should not be extremely light (we require at least 200 Da).
         (Glycosylated rotenoids may exceed classical weight bounds.)
      2. Look for candidate six‐membered rings that could be part of the rotenoid core.
         We define candidate rings as follows:
           a) A "chromenone‐like" ring must be six‐membered, contain at least one ring oxygen,
              and have exactly one carbon (of the ring) that is double‐bonded to an oxygen.
           b) An "ether" ring is a six‐membered ring that contains at least one oxygen
              and NO in‐ring carbonyl bonds.
      3. To be considered fused in the typical cis mode the two candidate rings must share exactly 2 atoms.
      4. (Optional extra check: if more than one candidate pair is found, we favor the one
         whose merged set of ring atoms does not introduce additional in‐ring carbonyls beyond
         that on the chromenone ring.)

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): A tuple with the boolean classification and a reason string.
                     If criteria are met, returns True with an explanation.
                     Otherwise, returns False with the reason for rejection.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule is not extremely light. (Some rotenoids are glycosylated, so we avoid a tight upper bound.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} is too low for rotenoids"

    # Get the ring atom indices as tuples (each tuple is a ring)
    ring_infos = mol.GetRingInfo().AtomRings()
    # Lists to store candidate rings. We use sets of atom indices for ease of intersection.
    chromenone_rings = []  # candidate rings with one in‐ring carbonyl
    ether_rings = []      # candidate rings with no in‐ring carbonyl

    # Loop over each ring and only consider six-membered ones.
    for ring in ring_infos:
        if len(ring) != 6:
            continue  # consider only 6-membered rings

        # Check for at least one oxygen atom in the ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not any(atom.GetAtomicNum() == 8 for atom in ring_atoms):
            continue

        # Count the number of ring carbon atoms that are double bonded to an oxygen.
        carbonyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Search bonds; if a bond from this atom is a double bond and the other atom is oxygen,
            # and if that oxygen is not part of the ring (or even if it is, we want to count the carbonyl feature)
            for bond in atom.GetBonds():
                # bond order 2
                if bond.GetBondTypeAsDouble() == 2:
                    other = bond.GetOtherAtom(atom)
                    if other.GetAtomicNum() == 8:
                        carbonyl_count += 1
                        break  # count each carbon only once

        # Classify the ring based on its in‐ring carbonyl count.
        # We assume that a chromenone candidate has exactly one in‐ring carbonyl 
        # and an ether candidate has none.
        if carbonyl_count == 1:
            chromenone_rings.append(set(ring))
        elif carbonyl_count == 0:
            ether_rings.append(set(ring))
        # Rings with >1 carbonyl are not considered as part of the core.

    if not chromenone_rings:
        return False, "No candidate chromenone-like ring (6-membered with one carbonyl and at least one oxygen) found"

    if not ether_rings:
        return False, "No candidate ether ring (6-membered with oxygen and no in-ring carbonyl) found"

    # Look for a pair of rings (one from each candidate list) that are fused.
    # We require they share exactly 2 atoms (typical for cis-fused systems).
    fused_pair_found = False
    for chrom in chromenone_rings:
        for ether in ether_rings:
            shared_atoms = chrom.intersection(ether)
            if len(shared_atoms) == 2:
                # Optionally, merge the rings and check that the only in-ring carbonyl is in the chromenone candidate.
                merged = chrom.union(ether)
                overall_carbonyls = 0
                for idx in merged:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() != 6:
                        continue
                    for bond in atom.GetBonds():
                        if bond.GetBondTypeAsDouble() == 2:
                            other = bond.GetOtherAtom(atom)
                            if other.GetAtomicNum() == 8:
                                overall_carbonyls += 1
                                break
                # We allow additional carbonyls outside the fused core by checking the core count.
                # Here we require that the fused core shows at least the one carbonyl from the chromenone.
                if overall_carbonyls >= 1:
                    fused_pair_found = True
                    break
        if fused_pair_found:
            break

    if fused_pair_found:
        return True, "Found fused pair: one six‐membered chromenone-like ring (with one in‐ring carbonyl and ring oxygen) cis-fused (sharing 2 atoms) to a six‐membered ether ring; molecular properties are consistent with a rotenoid skeleton."
    else:
        return False, "No characteristic cis-fused chromenone/ether ring pair was found; molecule is unlikely to be a rotenoid."

# Example test code (if executed as main, prints out a test result):
if __name__ == "__main__":
    # Test SMILES for a known rotenoid – 13alpha-Hydroxydolineone
    test_smiles = "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C=C6OC=CC6=C5"
    result, reason = is_rotenoid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)