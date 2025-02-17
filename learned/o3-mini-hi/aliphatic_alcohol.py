"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: An aliphatic alcohol, defined as 'An alcohol derived from an aliphatic compound.'

Improvement rationale:
- Instead of disqualifying a candidate –OH if any neighbor is aromatic, we now only consider
  candidates where the –OH is attached to a saturated (sp3) carbon that is non‐aromatic and is not
  part of any ring. (A cyclic –OH, such as in sugars or steroids, can often be part of a more complex 
  structure that should not be classified as a simple aliphatic alcohol.)
- We then perform a search over the molecule’s carbon network (restricted to non‐aromatic, sp3 carbons
  not in rings) to check if there is a long enough aliphatic chain (we use a threshold of at least 6 carbons)
  that “supports” the classification.
- Only if a candidate –OH is found on such an aliphatic chain, we classify the molecule as an aliphatic alcohol.
  
This approach helped us avoid false positives where a lone –OH on a (possibly benzylic) carbon triggers a classification,
and it also recovers some of the false negatives.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    Here, aliphatic alcohols are defined as having at least one hydroxyl (-OH) group attached to a saturated 
    (sp3), non-aromatic carbon that is not in a ring. Furthermore, that carbon must be part of a contiguous aliphatic 
    chain (of at least 6 carbons) to ensure that the alcohol is derived from an essentially aliphatic compound.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Helper function: perform DFS from a starting carbon atom over restricted aliphatic (sp3, non-aromatic, non-ring) carbons.
    def dfs_aliphatic_chain(atom, visited):
        visited.add(atom.GetIdx())
        length = 1  # count this atom
        for nbr in atom.GetNeighbors():
            # Only follow carbons that are saturated (sp3), non-aromatic, not in a ring, and not visited already.
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                if nbr.GetIdx() not in visited:
                    branch_length = dfs_aliphatic_chain(nbr, visited)
                    # We take maximum path length among branches
                    length = max(length, 1 + branch_length)
        return length

    # Iterate over oxygen atoms (possible hydroxyl groups)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Check if it is part of an -OH group (has at least one H either explicitly or implicitly)
        if atom.GetTotalNumHs() < 1:
            continue  # not an -OH group
        # Look at neighbors: we need one of them to be a carbon candidate
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                # Focus on the candidate carbon: it must be sp3, non-aromatic and not part of a ring.
                if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or nbr.GetIsAromatic() or nbr.IsInRing():
                    continue
                # Now, from this candidate carbon, search for a long contiguous chain of aliphatic carbons.
                chain_length = dfs_aliphatic_chain(nbr, set())
                if chain_length >= 6:
                    return True, ("Found -OH group attached to a non-ring, saturated (sp3, non-aromatic) carbon "
                                  "that is part of a contiguous aliphatic chain (chain length = {}).".format(chain_length))
    return False, ("No qualifying -OH group found on a non-ring, sp3 carbon that is part of a sufficiently long "
                   "aliphatic chain (>= 6 carbons).")
    
# (Optional) Testing examples – these are not required in the final program.
if __name__ == "__main__":
    test_smiles = [
        "O=C1OC([C@@H](O)\\C=C/C=C/C)CC1",  # Sapinofuranone A (True)
        "CCCCCCC(C)O",                     # octan-2-ol (True)
        "O=C(OC)/C=C/CC(O)CCCCCC(O)C",       # Cladosporester A (True)
        "OCCCCCC/C=C/CCCC",                # 7E-dodecen-1-ol (True)
        "O[C@H](C1=CC=C(CO)C=C1)CCCCC",     # 1-[4-(Hydroxymethyl)phenyl]hexan-1-ol (should be True)
        "OC1=CC=CC=C1",                    # Phenol: -OH on aromatic carbon (False)
    ]
    
    for sm in test_smiles:
        is_valid, reason = is_aliphatic_alcohol(sm)
        print("SMILES:", sm)
        print("Classification:", is_valid)
        print("Reason:", reason)
        print("-" * 60)