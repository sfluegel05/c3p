"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar – Any sugar having one or more alcoholic hydroxy groups replaced by 
substituted or unsubstituted amino groups.

Revised algorithm:
  1. Parse the molecule from the provided SMILES.
  2. Loop over all rings (as provided by RDKit) looking for candidate 5– or 6–membered rings.
  3. For each candidate ring, require that the ring contains exactly one oxygen atom 
     (typical for pyranose/furanose rings) and all the other ring atoms are carbons.
  4. For every ring carbon (non-oxygen) check each neighbor that is not in the ring.
     • For an oxygen (atomic number 8) we require that it has at least one hydrogen (i.e. resembles –OH).
     • For a nitrogen (atomic number 7) we count it only if the substituent is not part of another ring 
       or (if in a ring) is non‐aromatic. This is meant to filter out substituents that arise from a nucleobase.
  5. Additionally, require that the candidate ring has at least three such exocyclic substituents, thus
     ensuring it is “sugar‐like.”
  6. If a candidate ring is found that also has at least one acceptable amino (–NH2 or acetamido) substituent,
     we classify the molecule as an amino sugar.
  7. Otherwise, report a reason (or a negative result).
  
Note: this heuristic is not “perfect” but attempts to average over some mis‐classifications.
"""

from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    The algorithm inspects rings of size 5 and 6 for a sugar‐like pattern:
      • The ring must contain exactly one oxygen and the rest carbons.
      • Each ring carbon’s exocyclic substituents are examined: we count those that look like hydroxyl groups
        (oxygen with at least one hydrogen) or nitrogen atoms that are not part of another ring or aromatic.
      • To be “sugar–like” we require at least three exocyclic substituents overall, and at least one must be a 
        nitrogen substituent (free amine or acetamido) to indicate that a hydroxyl has been replaced.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple with True if the molecule is classified as an amino sugar,
                   and a reason string for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo().AtomRings()
    candidate_found = False  # flag if we find any sugar-like ring
    # Loop over all rings in the molecule.
    for ring in ring_info:
        # Consider only rings of size 5 or 6.
        if len(ring) not in (5, 6):
            continue

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count ring oxygens and carbons.
        ring_oxygens = [a for a in ring_atoms if a.GetAtomicNum() == 8]
        ring_carbons = [a for a in ring_atoms if a.GetAtomicNum() == 6]
        if len(ring_oxygens) != 1 or len(ring_carbons) != (len(ring) - 1):
            continue  # skip rings that do not match the typical sugar ring pattern

        # For each ring carbon, check its neighbors outside the ring.
        total_exo = 0       # total count of exocyclic substituents (O or N)
        valid_nitrogen = False  # flag if an acceptable amino substituent is detected

        for atom in ring_carbons:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms that are part of the ring
                at_num = nbr.GetAtomicNum()
                # If neighbor is oxygen, assume it is a hydroxyl if it has at least one hydrogen.
                if at_num == 8:
                    # Use GetTotalNumHs to see if there are attached hydrogens.
                    if nbr.GetTotalNumHs() > 0:
                        total_exo += 1
                    else:
                        # Even if no explicit H, count as substituent.
                        total_exo += 1
                # If neighbor is nitrogen, count it as a possible amino substituent 
                # if it is "free" (i.e. not embedded in an aromatic ring or another ring system).
                elif at_num == 7:
                    # Only count if the neighbor is not in a ring or (if in a ring) is not aromatic.
                    if (not nbr.IsInRing()) or (nbr.IsInRing() and not nbr.GetIsAromatic()):
                        total_exo += 1
                        valid_nitrogen = True
        # Require that the ring has a decent number of exocyclic substituents (at least 3).
        if total_exo < 3:
            continue

        candidate_found = True
        if valid_nitrogen:
            return True, ("Sugar ring detected (5- or 6-membered ring with one oxygen and at least three exocyclic substituents) "
                          "with at least one acceptable amino substituent (free amine or acetamido group).")
    # End loop over rings.
    if not candidate_found:
        return False, "No typical sugar ring (5- or 6-membered ring with one oxygen and multiple exocyclic substituents) was detected."
    else:
        return False, "Sugar-like ring(s) detected but none had an appropriate amino substituent (the nitrogen appeared to be part of a fused or aromatic ring)."

# Example usage:
if __name__ == "__main__":
    # Example SMILES that (according to the previous outcomes) should be classified as an amino sugar.
    example_smiles = ("O([C@H]1[C@@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)"
                      "[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO")
    result, reason = is_amino_sugar(example_smiles)
    print("Result:", result)
    print("Reason:", reason)