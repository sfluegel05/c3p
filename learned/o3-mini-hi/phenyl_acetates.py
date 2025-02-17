"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
Definition: An acetate ester obtained by the formal condensation of the carboxy group of acetic acid 
with the hydroxy group of any phenol.
A valid phenyl acetate must contain an acyclic acetoxy group –O–C(=O)CH3 that is directly attached
to a benzene ring at the site where a simple (phenolic) substituent is present.
We first look for the SMARTS substructure "[c:1][O:2][C:3](=O)[CH3]".
Then for each match we:
  - Ensure the ester O–C(=O) bond is not in a ring (i.e. not a lactone).
  - Identify a benzene ring (6-membered, all carbons aromatic) that contains the aromatic carbon [c:1].
  - Check that the substituents (other than the acetoxy group) attached to that aromatic carbon are “small”
    (have a heavy-atom count ≤ 5).
If a candidate match passes these tests, we classify the molecule as a phenyl acetate.
"""

from rdkit import Chem

def _substituent_heavy_atom_count(mol, start_idx, exclude):
    """
    Starting from the neighbor atom (start_idx), count the number of heavy atoms (C, N, O, etc)
    in the substituent branch. We perform a DFS not crossing any atom index in the "exclude" set.
    (Hydrogens are not counted.)
    """
    count = 0
    stack = [start_idx]
    visited = set()
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        atom = mol.GetAtomWithIdx(current)
        # Count heavy atoms (atomic number > 1)
        if atom.GetAtomicNum() > 1:
            count += 1
        # Get neighbors that are not in the exclude list and not visited.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in visited and nbr_idx not in exclude:
                stack.append(nbr_idx)
    return count

def is_phenyl_acetates(smiles: str):
    """
    Determines whether a molecule is a phenyl acetate as defined (an acetate ester
    derived from a phenol) based on its SMILES string.
    
    This function first looks for the acetoxy substructure: an aromatic carbon bonded to an oxygen
    that is in turn bonded to a C(=O)CH3 moiety. It further ensures that:
      1. The ester bond (O–C(=O)) is not in a ring (ruling out lactones).
      2. The aromatic carbon (that would be the original –OH position) belongs to a benzene ring
         (i.e. a 6-membered ring where every atom is carbon and marked aromatic).
      3. That aromatic carbon does not carry a “large” substituent (aside from the acetoxy group)
         – we compute (roughly) the heavy-atom count of any non–ring branch attached to that carbon.
         (The assumption is that in a true phenyl acetate the site of acetylation comes from a simple phenol.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a phenyl acetate, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for an acetoxy group attached to an aromatic carbon.
    acetate_pattern = Chem.MolFromSmarts("[c:1][O:2][C:3](=O)[CH3]")
    if acetate_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No phenyl acetate substructure ([c][O][C](=O)C) found in the molecule."
    
    # Get ring info for later ring membership tests
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Process each match candidate
    for match in matches:
        # match structure: (aromatic carbon, oxygen, carbonyl carbon, methyl carbon)
        aro_idx, o_idx, carbonyl_idx, methyl_idx = match
        
        # Check that the ester bond between the oxygen (o_idx) and carbonyl (carbonyl_idx) is NOT in a ring.
        in_ring = False
        for ring in ring_info:
            if o_idx in ring and carbonyl_idx in ring:
                in_ring = True
                break
        if in_ring:
            # This acetoxy group is in a ring (lactone) so skip it.
            continue
        
        # Verify that the aromatic carbon (aro_idx) is indeed aromatic.
        aro_atom = mol.GetAtomWithIdx(aro_idx)
        if not aro_atom.GetIsAromatic():
            continue
        
        # Try to find a 6-membered ring (benzene ring) that contains aro_idx.
        benzene_ring = None
        for ring in ring_info:
            if aro_idx in ring and len(ring) == 6:
                # Check that each atom in the ring is a carbon and is aromatic.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() 
                       for idx in ring):
                    benzene_ring = set(ring)
                    break
        if benzene_ring is None:
            # The aromatic atom is not in a simple benzene ring.
            continue
        
        # For the candidate benzene ring, check the substituents on the aromatic carbon (aro_idx)
        # apart from the acetoxy oxygen (o_idx) and the ring connections.
        aro_neighbors = [nbr.GetIdx() for nbr in aro_atom.GetNeighbors()]
        # Build a set of indices to exclude from the substituent search: the entire ring.
        exclude = set(benzene_ring)
        # Also exclude the acetoxy oxygen.
        exclude.add(o_idx)
        # Now, for each neighbor of the aromatic carbon that is not in the benzene ring
        # (i.e. substituents off the ring at that position) compute its heavy atom count.
        large_substituent_found = False
        for nbr_idx in aro_neighbors:
            if nbr_idx in benzene_ring or nbr_idx == o_idx:
                continue
            count = _substituent_heavy_atom_count(mol, nbr_idx, exclude)
            # If any substituent branch has more than 5 heavy atoms, consider it too big.
            if count > 5:
                large_substituent_found = True
                break
        
        if large_substituent_found:
            # This match seems to be on a benzene ring that is substituted with a large fragment;
            # likely not a simple phenol-derived acetate.
            continue
        
        # If we reach here, the candidate acetoxy group is acyclic and attached to a benzene ring
        # with only small substituents at the attachment point.
        return True, "Molecule contains an acyclic phenyl acetate moiety (acetylated phenolic oxygen) with acceptable substituents."
    
    # If no candidate match passes all tests
    return False, "Acetoxy group found but did not meet phenyl acetate criteria (either in a ring or attached to an over‐substituted aromatic center)."

# Example usage:
if __name__ == "__main__":
    # Test examples (one true positive and one false positive from the outcomes)
    test_smiles_true = "CC(=O)Oc1ccccc1"  # phenyl acetate
    result, reason = is_phenyl_acetates(test_smiles_true)
    print("Test (phenyl acetate):", result, reason)
    
    test_smiles_false = "CCOC(=O)C1=C(N=C2N(C1C3=CC(=C(C=C3)OC(=O)C)OC)C(=O)C(S2)C)C"
    result, reason = is_phenyl_acetates(test_smiles_false)
    print("Test (false positive candidate):", result, reason)