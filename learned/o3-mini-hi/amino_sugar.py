"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar – Any sugar having one or more alcoholic hydroxy groups
replaced by substituted or unsubstituted amino groups.
Improved algorithm:
  1. Parse the molecule.
  2. Search for candidate rings of size 5 or 6 that have exactly one ring oxygen.
  3. For each candidate ring, check that the “sugar‐like” pattern is met by requiring
     that many of the ring carbons (non‐oxygen atoms) have exocyclic substituents that are either oxygen or nitrogen.
  4. Then check if at least one such substituent is a nitrogen (allowing for –NH2 or acetamido groups),
     which we interpret as a hydroxyl replacement.
  5. Otherwise, return a negative result with an appropriate message.
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string. 
    The algorithm looks for a ring (5 or 6 membered) containing exactly one oxygen
    (as is typical in furanoses or pyranoses) and expects that several of the ring carbons 
    carry exocyclic oxygen substituents (the hydroxyl groups found in sugars). 
    If one or more of these hydroxyl substituents is replaced by nitrogen (or an acylated nitrogen),
    we classify the structure as an amino sugar.
    
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
    candidate_found = False

    # Iterate over each ring in the molecule.
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        
        # Count the oxygen atoms within the ring.
        ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygens) != 1:
            continue  # typical sugar ring should have one oxygen in the ring
        
        # Now, consider the ring atoms that are carbon.
        ring_carbons = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        if not ring_carbons:
            continue

        # For a typical sugar, several ring carbons have exocyclic substituents that are O (or N).
        subs_count = 0
        amino_found = False
        
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            # Count substituents not in the ring.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                # Accept oxygen or nitrogen substituents as expected in normal sugars.
                if neighbor.GetAtomicNum() in (7, 8):
                    subs_count += 1
                    # For amino substitution, the substituent nitrogen should be directly attached.
                    if neighbor.GetAtomicNum() == 7:
                        amino_found = True
        
        # Require that the ring shows enough “sugar-like” substitution.
        # Typical sugars have multiple -OH groups so we expect at least 2 exocyclic O or N groups.
        if subs_count < 2:
            continue
        
        candidate_found = True  # we have a sugar-like ring candidate
        
        if amino_found:
            return True, "Sugar ring detected with at least one alcoholic hydroxyl replaced by an amino (or acylated amino) group."

    if not candidate_found:
        return False, "No typical sugar ring (5- or 6-membered ring with one oxygen and multiple exocyclic substituents) was detected."
    else:
        return False, "Sugar-like ring detected but no amino substitution on a ring carbon was found."

# Example usage:
if __name__ == "__main__":
    example_smiles = "O([C@H]1[C@@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO"
    result, reason = is_amino_sugar(example_smiles)
    print("Result:", result)
    print("Reason:", reason)