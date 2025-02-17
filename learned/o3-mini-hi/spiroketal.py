"""
Classifies: CHEBI:72600 spiroketal
"""
#!/usr/bin/env python
"""
Classifies: spiroketal
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.

A true spiroketal features a sp³ carbon that sits at the junction of two rings and that
only bridges the two rings. In a cyclic ketal the spiro carbon is tetravalent, with exactly
two oxygen substituents and two carbon substituents. Furthermore, in a true spiroketal each of
the two rings (that share only the spiro center) should include one (and only one)
oxygen substituent. This program checks those criteria.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.

    A candidate spiroketal center must:
      - Be a carbon (atomic number 6) that is sp3 hybridized.
      - Be tetravalent (exactly 4 bonds) with exactly 2 oxygen and 2 carbon substituents.
      - Belong to at least two rings.
      - Have at least one pair of rings that share only the candidate atom.
      - Additionally, for that pair of rings, one oxygen neighbor should appear in each ring;
        that is, the two oxygen substituents must “separate” between the two rings.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Explanation for the classification.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information (a tuple of atom index tuples)
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings detected in the molecule"

    # Loop over all atoms searching for candidate spiroketal carbon
    for atom in mol.GetAtoms():
        # Consider only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue

        # Check that the atom is sp3-hybridized.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue

        atom_idx = atom.GetIdx()

        # Candidate must appear in at least 2 rings.
        atom_rings = [set(r) for r in rings if atom_idx in r]
        if len(atom_rings) < 2:
            continue

        # Check for a pair of rings that share only this candidate.
        candidate_ring_pair_found = False
        candidate_ring_pair = None
        for i in range(len(atom_rings)):
            for j in range(i+1, len(atom_rings)):
                if atom_rings[i].intersection(atom_rings[j]) == {atom_idx}:
                    candidate_ring_pair_found = True
                    candidate_ring_pair = (atom_rings[i], atom_rings[j])
                    break
            if candidate_ring_pair_found:
                break
        if not candidate_ring_pair_found:
            continue

        # Now check the ketal pattern:
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 4:
            continue  # Must be tetravalent

        oxygen_count = 0
        carbon_count = 0
        # Record neighbor indices by element type.
        oxygen_neighbors = []
        carbon_neighbors = []
        for nbr in neighbors:
            anum = nbr.GetAtomicNum()
            if anum == 8:
                oxygen_count += 1
                oxygen_neighbors.append(nbr.GetIdx())
            elif anum == 6:
                carbon_count += 1
                carbon_neighbors.append(nbr.GetIdx())
            # Other atoms are allowed but we focus on oxygen vs carbon here.
        if oxygen_count != 2 or carbon_count != 2:
            continue

        # For the candidate ring pair, check that each ring contains one, and only one,
        # distinct oxygen neighbor.
        ring1, ring2 = candidate_ring_pair
        ox_in_ring1 = [ox for ox in oxygen_neighbors if ox in ring1]
        ox_in_ring2 = [ox for ox in oxygen_neighbors if ox in ring2]

        # We do not want the same oxygen to be in both rings.
        if not ox_in_ring1 or not ox_in_ring2:
            continue
        if set(ox_in_ring1).intersection(set(ox_in_ring2)):
            # If an oxygen appears in both rings, the rings share more than the candidate.
            continue

        # If the criteria are satisfied, then we classify the molecule as spiroketal.
        return True, (f"Found spiroketal center at atom index {atom_idx} with "
                      f"{oxygen_count} oxygen and {carbon_count} carbon substituents, "
                      f"and the oxygen atoms distribute uniquely in two rings "
                      f"(ring pair with only common atom {atom_idx}).")

    return False, "No spiroketal pattern found in the molecule"

# Example usage:
if __name__ == "__main__":
    # Test with one of the examples (nocamycin E)
    test_smiles = "[C@]12([C@@H]([C@](O[C@]3(O1)[C@H]([C@H](O[C@]3(CC2=O)C)C)C(=O)O)([C@H](C)/C=C(\\C)/C=C/C(=C/4\\C(NCC4=O)=O)/O)[H])C)[H]"
    result, reason = is_spiroketal(test_smiles)
    print(result, reason)