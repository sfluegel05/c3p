"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide – Any aromatic amide obtained by acylation of aniline.

An anilide is defined as an aromatic amide where the nitrogen (originally from aniline) 
is connected directly to a benzene ring and to an acyl group. In other words, the structure 
should include a benzene ring attached to a nitrogen that is bonded to a carbonyl group.
  
This version uses two SMARTS patterns (with a more permissive ‘N’ instead of “NH”) together 
with an extra check on the candidate nitrogen:
  1. Pattern1: "c1ccccc1NC(=O)" – a benzene ring directly attached to an N that is acylated.
  2. Pattern2: "NC(=O)c1ccccc1" – the acyl group comes before the benzene ring.
  
For each match the candidate N atom is validated by checking that:
  • It is not part of any ring.
  • It has exactly 2 heavy-atom (non-hydrogen) neighbors.
  • One neighbor is an aromatic carbon (from a benzene ring) 
    and the other is a carbon that is part of a carbonyl (i.e. has a double bond to oxygen).
If a match passes these tests the molecule is classified as an anilide.
"""

from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.

    An anilide is defined as an aromatic amide originating from aniline (i.e. the N from aniline 
    is acylated). The key requirement is that the amide nitrogen must be directly bonded to an aromatic 
    benzene ring carbon and to a carbonyl carbon. This function uses two SMARTS patterns with extra
    checks on the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an anilide, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns. We do not use an explicit hydrogen on N, so that tertiary N's are accepted.
    pattern1 = Chem.MolFromSmarts("c1ccccc1NC(=O)")
    pattern2 = Chem.MolFromSmarts("NC(=O)c1ccccc1")
    
    # Helper function to validate that a candidate nitrogen is likely derived from aniline.
    def valid_anilide_nitrogen(n_atom):
        # The aniline N should not be part of a ring.
        if n_atom.IsInRing():
            return False
        
        # Count heavy neighbors (non-hydrogen)
        heavy_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 2:
            return False
        
        aromatic_ok = False  # one neighbor should be an aromatic carbon (from benzene)
        acyl_ok = False      # the other neighbor should be a carbonyl carbon (showing a double bond to O)
        
        for nbr in n_atom.GetNeighbors():
            # Check for aromatic neighbor: carbon and aromatic flag True.
            if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                aromatic_ok = True
            # Check for acyl carbon: carbon having a double-bond to oxygen.
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            acyl_ok = True
        return aromatic_ok and acyl_ok

    # Check using pattern1. In the SMARTS "c1ccccc1NC(=O)", the N is expected at the 7th position.
    matches1 = mol.GetSubstructMatches(pattern1)
    for match in matches1:
        # Expect at least 7 atoms in the match (indices 0-6 for the N to be at 6)
        if len(match) < 7:
            continue
        n_idx = match[6]
        n_atom = mol.GetAtomWithIdx(n_idx)
        if valid_anilide_nitrogen(n_atom):
            return True, "Found anilide substructure: benzene ring attached to N acyl group (pattern1)"
    
    # Check using pattern2. In the SMARTS "NC(=O)c1ccccc1", the N should be the first atom (index 0).
    matches2 = mol.GetSubstructMatches(pattern2)
    for match in matches2:
        if len(match) < 1:
            continue
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        if valid_anilide_nitrogen(n_atom):
            return True, "Found anilide substructure: N acyl group linked to benzene (pattern2)"
    
    return False, "No valid anilide substructure found"

# Example usage. You can run a few test cases.
if __name__ == "__main__":
    # Example: one of the true positives
    test_smiles = "CC(=O)NC1=CC=C(C=C1)C(=S)NCC2=CC=CO2"
    result, reason = is_anilide(test_smiles)
    print("Is anilide:", result)
    print("Reason:", reason)