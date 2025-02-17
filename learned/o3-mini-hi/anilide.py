"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide – Any aromatic amide obtained by acylation of aniline.
An anilide is defined as an aromatic amide where the amide nitrogen is derived from aniline,
i.e. the N comes from an NH2 group attached directly to a benzene ring, and subsequently acylated.
This version uses two SMARTS patterns plus additional checks on the candidate N:
  1. "c1ccccc1[NH]C(=O)" – the benzene ring is attached via –NH– to an acyl group.
  2. "NC(=O)c1ccccc1" – the acyl group comes first, but still the same N joins to a benzene ring.
For each match the candidate N atom is checked to ensure that it:
  • is not itself part of any ring (so it comes from aniline),
  • and has exactly 2 heavy-atom neighbors (one from the acyl group and one aromatic carbon).
If one of the matches meets the criteria the molecule is classified as an anilide.
"""

from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    
    An anilide is defined as an aromatic amide derived from the acylation of aniline.
    That is, the N in NC(=O)– must be directly attached to a benzene ring
    (i.e. a six-membered purely carbon aromatic system) and not be part of any ring itself.
    
    The function uses two SMARTS patterns:
      Pattern 1 ("c1ccccc1[NH]C(=O)"): benzene ring with an exocyclic NH that is acylated.
      Pattern 2 ("NC(=O)c1ccccc1"): the acyl portion comes first, but still the same connectivity.
    
    In addition, for each substructure match the following are verified:
      - The nitrogen (N) is not part of any ring.
      - The nitrogen is a secondary amine (has exactly 2 heavy atom, i.e. non-hydrogen, neighbors).
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as an anilide, False otherwise.
        str: Reason for classification.
    """
    
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns.
    # Pattern 1: The benzene ring is attached to an NH that is acylated.
    # In this SMARTS the match order is: 6 atoms for the benzene {indices 0-5}, then the N (index 6), then the carbonyl carbon (index 7)
    pattern1 = Chem.MolFromSmarts("c1ccccc1[NH]C(=O)")
    # Pattern 2: A slightly different direction; here the N (index 0) is acylated and then attached to a benzene ring.
    # The matching order will be: N (index 0), carbonyl C (index 1), O (index 2) and then benzene (indices 3-8)
    pattern2 = Chem.MolFromSmarts("NC(=O)c1ccccc1")
    
    # Helper function to check if a candidate N atom is “good”.
    def valid_nitrogen(n_atom):
        # The nitrogen should not be part of any ring,
        # and should have exactly 2 non-hydrogen neighbors.
        if n_atom.IsInRing():
            return False
        heavy_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        return len(heavy_neighbors) == 2
    
    # Try pattern 1: in this match, we expect the N to be at position 6.
    matches1 = mol.GetSubstructMatches(pattern1)
    for match in matches1:
        # Ensure match is long enough (should be 8 atoms).
        if len(match) < 7:
            continue
        # For pattern1, the N is at index 6.
        n_idx = match[6]
        n_atom = mol.GetAtomWithIdx(n_idx)
        if valid_nitrogen(n_atom):
            return True, "Found anilide substructure: benzene ring attached to NH acyl group (pattern1)"
    
    # Try pattern 2: here the nitrogen should occur at index 0.
    matches2 = mol.GetSubstructMatches(pattern2)
    for match in matches2:
        # For pattern2, expect the match to have at least 1 atom for the N.
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        if valid_nitrogen(n_atom):
            return True, "Found anilide substructure: acylated NH linked to benzene (pattern2)"
    
    return False, "No valid anilide substructure found"

# Example usage
if __name__ == "__main__":
    # Use one example from the true positives.
    test_smiles = "CC(=O)NC1=CC=C(C=C1)C(=S)NCC2=CC=CO2"
    result, reason = is_anilide(test_smiles)
    print("Is anilide:", result)
    print("Reason:", reason)