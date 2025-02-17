"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: Disaccharide
A disaccharide is defined as a compound in which two monosaccharides are joined by a glycosidic bond.
This improved algorithm works in two stages:
  1) Identify candidate monosaccharide rings by matching SMARTS for pyranose (6-membered with one oxygen)
     and furanose (5-membered with one oxygen) rings.
  2) Identify a glycosidic linkage by finding an oxygen atom (not in either candidate ring) that is bonded
     to a carbon in one ring and a carbon in the other.
If both conditions are met, the molecule is classified as a disaccharide.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    
    A disaccharide should contain exactly two monosaccharide rings (typically pyranose or furanose)
    connected by a glycosidic linkage. Here we:
      1) Use SMARTS patterns to identify candidate sugar rings.
      2) Look for an oxygen linker (not part of the rings) that bonds to carbons in both rings.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a disaccharide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    candidate_rings = []
    
    # Define SMARTS patterns for pyranose and furanose rings:
    # Pyranose: 6-membered ring with 5 carbons and 1 oxygen.
    pyranose_smarts = "[C;R6][C;R6][C;R6][C;R6][C;R6][O;R6]"
    # Furanose: 5-membered ring with 4 carbons and 1 oxygen.
    furanose_smarts = "[C;R5][C;R5][C;R5][C;R5][O;R5]"
    
    pyranose = Chem.MolFromSmarts(pyranose_smarts)
    furanose = Chem.MolFromSmarts(furanose_smarts)
    
    # Get matches for each pattern. Each match is a tuple of atom indices.
    pyranose_matches = mol.GetSubstructMatches(pyranose)
    furanose_matches = mol.GetSubstructMatches(furanose)
    
    # Collect candidate rings as sets of atom indices.
    for match in pyranose_matches:
        candidate_rings.append(set(match))
    for match in furanose_matches:
        candidate_rings.append(set(match))
        
    # Remove duplicate ring matches
    unique_rings = []
    seen = set()
    for ring in candidate_rings:
        fs = frozenset(ring)
        if fs not in seen:
            seen.add(fs)
            unique_rings.append(ring)
            
    # We expect exactly 2 monosaccharide rings for a disaccharide.
    if len(unique_rings) != 2:
        return False, f"Found {len(unique_rings)} candidate sugar ring(s); exactly 2 are needed for a disaccharide."
    
    ring1, ring2 = unique_rings
    
    # To find a glycosidic linkage, search for an oxygen atom (outside the candidate rings)
    # that is directly connected to at least one carbon atom from ring1 and one carbon atom from ring2.
    glyco_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # only consider oxygen atoms
        # Skip oxygens already belonging to one of the sugar rings (they are part of the ring template)
        if atom.GetIdx() in ring1 or atom.GetIdx() in ring2:
            continue
        # Check neighbors: typically the bridging oxygen connects to carbons in each sugar.
        neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # If it connects to at least one carbon in ring1 and at least one in ring2, regard it as a glycosidic oxygen.
        if any(nbr.GetIdx() in ring1 for nbr in neighbor_carbons) and any(nbr.GetIdx() in ring2 for nbr in neighbor_carbons):
            glyco_found = True
            break
    
    if glyco_found:
        return True, "Contains exactly two monosaccharide rings joined by a glycosidic bond."
    else:
        return False, "No glycosidic linkage bridging the two sugar rings was found."

# Example usage (for testing):
if __name__ == "__main__":
    # Provided test disaccharide: alpha-L-Fucp-(1->6)-alpha-D-Glcp
    test_smiles = "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)CO[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C"
    result, reason = is_disaccharide(test_smiles)
    print(result, reason)