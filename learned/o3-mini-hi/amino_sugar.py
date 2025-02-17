"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar – Any sugar having one or more alcoholic hydroxy groups
replaced by substituted or unsubstituted amino groups.
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    The algorithm does the following:
      1. Parse the SMILES string.
      2. Identify candidate sugar rings. Many sugars are cyclic ethers (furanoses (5-membered) 
         or pyranoses (6-membered)) that contain exactly one ring oxygen.
      3. For each such candidate ring, check the ring carbons for substituents.
         Under normal circumstances these carbons carry hydroxyl groups.
         If one or more of these substituents is a nitrogen (atomic number 7), then we consider
         that the corresponding -OH is replaced by an amino (or acetamido, etc.) group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule appears to be an amino sugar, False otherwise.
        str: A brief explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring info from molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_candidate_found = False
    amino_substitution_found = False

    # Iterate over each ring in the molecule.
    for ring in ring_info:
        # Many common sugar rings are five or six members
        if len(ring) not in (5, 6):
            continue

        # Count number of oxygen atoms in the ring.
        oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        # Typical furanose or pyranose rings have exactly one ring oxygen.
        if oxygen_count != 1:
            continue

        # We have a sugar-ring candidate
        sugar_candidate_found = True

        # For each atom index in the ring (usually the non-oxygen atoms are carbons)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We expect the ring to be carbons (the oxygen is part of the ring and not substituted)
            if atom.GetAtomicNum() != 6:
                continue
            # Check neighbor atoms that are not part of the ring.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                # If a substituent is a nitrogen, we interpret that as an amino replacement.
                if neighbor.GetAtomicNum() == 7:
                    amino_substitution_found = True
                    # We can return as soon as we find one substituted amino group.
                    return True, "Sugar ring detected with at least one alcoholic hydroxyl replaced by an amino group."

    if not sugar_candidate_found:
        return False, "No typical sugar ring (5- or 6-membered ring with one oxygen) was detected."

    if not amino_substitution_found:
        return False, "Sugar-like ring found but no amino substitution on a ring carbon was detected."

    return False, "Unknown reason"
    
# Example usage:
if __name__ == "__main__":
    # An example SMILES string taken from one of the provided examples.
    test_smiles = "O([C@H]1[C@@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO"
    result, reason = is_amino_sugar(test_smiles)
    print("Result:", result)
    print("Reason:", reason)