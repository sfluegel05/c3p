"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid
with the hydroxy group of any phenol.
A valid phenyl acetate must contain an acyclic acetoxy group –O–C(=O)CH3 directly attached to a benzene ring.
The strategy is:
  1. Look for the substructure "[c:1][O:2][C:3](=O)[CH3]".
  2. For each match, ensure that:
     - The ester bond (O–C(=O)) is not in a ring (to rule out lactones).
     - The aromatic carbon ([c:1]) belongs to a benzene ring (a 6-membered ring of aromatic carbons).
     - The substituents attached to that aromatic carbon (other than the acetoxy group) are small (heavy-atom count ≤ 5).
"""

from rdkit import Chem

def _substituent_heavy_atom_count(mol, start_idx, exclude):
    """
    Count the number of heavy atoms (non-hydrogen atoms) in the branch starting at start_idx,
    not crossing any atoms in the exclude set.
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
        if atom.GetAtomicNum() > 1:
            count += 1
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in visited and nbr_idx not in exclude:
                stack.append(nbr_idx)
    return count

def is_phenyl_acetates(smiles: str):
    """
    Determines whether a molecule is a phenyl acetate based on its SMILES string.
    
    The molecule must contain an acyclic acetoxy substructure where an acetoxy group (O–C(=O)CH3)
    is directly attached to an aromatic carbon of a benzene ring. In addition, the substituents 
    on the aromatic carbon must be small (heavy-atom count ≤ 5) so that it is consistent with 
    an acetylated simple phenol.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a phenyl acetate, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the acetoxy group attached to an aromatic carbon.
    acetate_pattern = Chem.MolFromSmarts("[c:1][O:2][C:3](=O)[CH3]")
    if acetate_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No phenyl acetate substructure ([c][O][C](=O)C) found."
    
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Process each candidate match.
    for match in matches:
        # Ensure the match has exactly 4 atoms. If not, skip this match.
        if len(match) != 4:
            continue
        aro_idx, o_idx, carbonyl_idx, methyl_idx = match
        
        # Check that the ester bond (O–C(=O)) is not in a ring (i.e. not a lactone).
        in_ring = False
        for ring in ring_info:
            if o_idx in ring and carbonyl_idx in ring:
                in_ring = True
                break
        if in_ring:
            continue  # Skip if the ester is part of a ring.
        
        # Verify that the aromatic carbon (aro_idx) is indeed aromatic.
        aro_atom = mol.GetAtomWithIdx(aro_idx)
        if not aro_atom.GetIsAromatic():
            continue
        
        # Find a benzene ring (6-member aromatic carbons) that contains this atom.
        benzene_ring = None
        for ring in ring_info:
            if aro_idx in ring and len(ring) == 6:
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic()
                       for idx in ring):
                    benzene_ring = set(ring)
                    break
        if benzene_ring is None:
            continue  # Aromatic atom not in a simple benzene ring.
        
        # Check substituents off the aromatic carbon (other than the acetoxy oxygen and ring members).
        aro_neighbors = [nbr.GetIdx() for nbr in aro_atom.GetNeighbors()]
        exclude = set(benzene_ring)
        exclude.add(o_idx)
        large_substituent_found = False
        for nbr_idx in aro_neighbors:
            if nbr_idx in benzene_ring or nbr_idx == o_idx:
                continue
            count = _substituent_heavy_atom_count(mol, nbr_idx, exclude)
            if count > 5:
                large_substituent_found = True
                break
        if large_substituent_found:
            continue
        
        # Passed all checks; this molecule is classified as a phenyl acetate.
        return True, "Molecule contains an acyclic phenyl acetate moiety (acetylated phenolic oxygen) with acceptable substituents."
    
    return False, "Acetoxy group found but did not meet phenyl acetate criteria (either in a ring or attached to an over‐substituted aromatic center)."

# Example usage:
if __name__ == "__main__":
    # Testing a known phenyl acetate
    test_smiles_true = "CC(=O)Oc1ccccc1"  # phenyl acetate
    result, reason = is_phenyl_acetates(test_smiles_true)
    print("Test (phenyl acetate):", result, reason)
    
    # Testing one of the provided examples that previously triggered an error.
    test_smiles_error = "COc1cc(ccc1OC(C)=O)C(C)=O"  # 4-acetyl-2-methoxyphenyl acetate
    result, reason = is_phenyl_acetates(test_smiles_error)
    print("Test (4-acetyl-2-methoxyphenyl acetate):", result, reason)