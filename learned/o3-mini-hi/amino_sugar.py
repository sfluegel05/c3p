"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar – Any sugar having one or more alcoholic hydroxy groups
replaced by substituted or unsubstituted amino groups.

Revised algorithm:
  1. Parse the molecule.
  2. Look for candidate 5- or 6-membered rings that have exactly one ring oxygen 
     (i.e. the typical pyranose or furanose ring) and the rest are carbons.
  3. For each ring-carbon (i.e. non-ring-oxygen) inspect its exocyclic neighbors.
     • Count a neighbor if it is oxygen and likely an –OH (has at least one hydrogen),
       or if it is a nitrogen atom (which we count as the amino substitution; note that nitrogen
       attached to a carbonyl is taken as acetamido).
  4. Require that the candidate sugar ring has (a) at least two exocyclic substituents (since 
     most sugars are polyhydroxylated) and (b) that at least one exocyclic substituent is a nitrogen.
  5. If such a candidate ring is found, classify the structure as an amino sugar.
  6. Otherwise, return a negative result with a reason.
"""

from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    The algorithm inspects rings of size 5 and 6 for a sugar‐like pattern:
    • Exactly one ring oxygen and the remaining ring atoms are carbons.
    • Multiple exocyclic substituents (as expected for hydroxyl groups) and at least
      one of these substituents is a nitrogen atom (as a free amine or acetamido group)
      implying a hydroxyl replacement.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as an amino sugar, False otherwise.
      str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_found = False  # have we seen any candidate sugar-like rings?
    
    # Iterate over each ring in the molecule.
    for ring in ring_info:
        # Consider only rings of size 5 or 6.
        if len(ring) not in (5, 6):
            continue
            
        # Check that the ring has exactly one oxygen atom and the rest are carbons.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_oxygens = [a for a in ring_atoms if a.GetAtomicNum() == 8]
        ring_carbons = [a for a in ring_atoms if a.GetAtomicNum() == 6]
        if len(ring_oxygens) != 1 or len(ring_carbons) != (len(ring) - 1):
            continue  # not a typical sugar ring
        
        # For each ring carbon, check exocyclic (outside of ring) substituents.
        total_subs = 0    # count how many exocyclic substituents (O or N) overall
        amino_found = False  # flag if at least one substituent is nitrogen
        hydroxyl_subs = 0  # count potential hydroxyl (O with attached H)
        
        for a in ring_carbons:
            # Loop over neighbors: skip those in the ring.
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                atm_num = nbr.GetAtomicNum()
                # If neighbor is oxygen, we assume it is a hydroxyl if it has any hydrogen
                if atm_num == 8:
                    # Check neighbor explicit hydrogens or using GetTotalNumHs
                    if nbr.GetTotalNumHs() > 0:
                        hydroxyl_subs += 1
                        total_subs += 1
                    else:
                        # Might be involved in an ester or ether; still count as substituent.
                        total_subs += 1
                # If neighbor is nitrogen, count directly as amino substituent.
                elif atm_num == 7:
                    total_subs += 1
                    amino_found = True
        # To be sugar-like we expect at least two exocyclic substituents.
        if total_subs < 2:
            continue
        
        candidate_found = True  # found a sugar-like ring candidate
        if amino_found:
            return True, ("Sugar ring detected (5- or 6-membered ring with one oxygen and multiple exocyclic substituents) "
                          "with at least one substituent being nitrogen (free amine or acetamido group).")
    # End of candidate ring search.
    if not candidate_found:
        return False, "No typical sugar ring (5- or 6-membered ring with one oxygen and multiple exocyclic substituents) was detected."
    else:
        return False, "Sugar-like ring detected but no amino substitution on a ring carbon was found."

# Example usage:
if __name__ == "__main__":
    # This example SMILES is taken from one of the reported true positives.
    example_smiles = "O([C@H]1[C@@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO"
    result, reason = is_amino_sugar(example_smiles)
    print("Result:", result)
    print("Reason:", reason)