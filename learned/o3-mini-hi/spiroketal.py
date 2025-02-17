"""
Classifies: CHEBI:72600 spiroketal
"""
#!/usr/bin/env python
"""
Classifies: spiroketal
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
A true spiroketal features a sp³ carbon that sits at the junction of two rings and that
only bridges the two rings. In a cyclic ketal the spiro carbon is tetravalent, with exactly
two oxygen and two carbon substituents. Furthermore, in a true spiroketal each of the two rings
(that share only the spiro center) should include exactly one oxygen substituent (from the spiro center).
The new approach below refines the pairing of oxygen substituents with the individual rings.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.

    The algorithm looks for a candidate spiro carbon that:
      - Is a carbon (atomic number 6) that is sp3 hybridized.
      - Has exactly 4 neighbors (i.e. is tetravalent) with exactly 2 oxygen and 2 carbon substituents.
      - Belongs to at least two rings.
      - For each of its two oxygen neighbors, there is at least one ring (of minimal size, here >=4 atoms)
        that contains both the candidate and that oxygen.
      - There exists a pair of such rings (one for each oxygen) whose intersection is only the candidate atom.
        In other words, the two rings “meet” only at the candidate spiro center (thus forming a spiro junction).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Retrieve ring information as a list of sets of atom indices.
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if not rings:
        return False, "No rings detected in the molecule"

    # Loop over all atoms looking for candidate spiro centers.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue

        atom_idx = atom.GetIdx()
        # Candidate must appear in at least 2 rings.
        candidate_rings = [r for r in rings if atom_idx in r and len(r) >= 4]  # ignore very small rings
        if len(candidate_rings) < 2:
            continue

        # Check degree and neighbor element counts.
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 4:
            continue  # Must be tetravalent.
        oxygen_neighbors = []
        carbon_neighbors = []
        for nbr in neighbors:
            anum = nbr.GetAtomicNum()
            if anum == 8:
                oxygen_neighbors.append(nbr.GetIdx())
            elif anum == 6:
                carbon_neighbors.append(nbr.GetIdx())
        if len(oxygen_neighbors) != 2 or len(carbon_neighbors) != 2:
            continue

        # For each oxygen neighbor, compile the list of rings (from candidate_rings)
        # that contain both the candidate and that oxygen.
        ox_rings = {}
        for ox_idx in oxygen_neighbors:
            ox_rings[ox_idx] = [r for r in candidate_rings if ox_idx in r]
            if not ox_rings[ox_idx]:
                # If an oxygen neighbor is not in any candidate ring, cannot form spiroketal.
                break
        else:  # Only proceed if both oxygens have a ring.
            # Now try to pair one ring from the first oxygen with one ring from the second oxygen.
            for r1 in ox_rings[oxygen_neighbors[0]]:
                for r2 in ox_rings[oxygen_neighbors[1]]:
                    # Check that the two rings share only the candidate carbon.
                    if r1.intersection(r2) == {atom_idx}:
                        # We have found the correct pairing.
                        return True, (f"Found spiroketal center at atom index {atom_idx} with two oxygens "
                                       f"({oxygen_neighbors[0]} and {oxygen_neighbors[1]}) and two carbons "
                                       f"({carbon_neighbors[0]} and {carbon_neighbors[1]}). "
                                       f"Ring one (size {len(r1)}) and ring two (size {len(r2)}) share only the candidate atom.")
    return False, "No spiroketal pattern found in the molecule"

# Example usage:
if __name__ == "__main__":
    # Try testing with one of the known examples, for instance chrysopine:
    test_smiles = "NC(=O)CC[C@@H]1NCC2(OC[C@@H](O)[C@@H](O)[C@@H]2O)OC1=O"
    result, reason = is_spiroketal(test_smiles)
    print(result, reason)