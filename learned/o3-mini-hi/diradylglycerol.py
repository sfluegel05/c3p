"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition:
  Any lipid that is glycerol bearing two substituent groups – either acyl, alkyl, or alk-1-enyl –
  at any two of the three possible positions.
  
This revised heuristic:
  1. Parses the SMILES and rejects invalid input.
  2. Immediately excludes molecules containing phosphorus.
  3. Searches for a glycerol-like three-carbon backbone. In “true” glycerol the backbone carbons 
     are non‐ring sp3 carbons arranged linearly (the two terminal carbons being CH2 and the central one CH).
  4. For each backbone carbon we require exactly one oxygen neighbor (not part of the backbone).
     Then, for each such oxygen we decide if it is still a free –OH (if it has any attached hydrogen)
     or has been substituted. In the latter case (no H on the oxygen), we perform a depth-first search
     from that oxygen (excluding the backbone) to count connected carbon atoms. If the branch contains
     at least 3 carbons we call that substituent “valid.”
  5. If exactly two backbone positions are substituted (and the third remains free) we return positive.
  
Note: This is still a heuristic approach. Some molecules “look like” glycerol --
either in a diradylglycerol or in a larger sugar or lipid – so edge cases remain.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    
    A diradylglycerol is defined as a glycerol backbone (a linear chain of three sp3 carbons,
    with each carbon bearing one oxygen substituent) in which exactly two of the oxygen groups 
    have been replaced by acyl/alkyl/alk-1-enyl chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diradylglycerol, otherwise False.
        str: Explanation for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Exclude molecules containing phosphorus (P, atomic number 15)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid, not a diradylglycerol"
    
    # Helper function: count the number of carbon atoms connected (via a DFS) starting from atom,
    # excluding atoms in exclude_idxs. Only heavy (atomic num > 1) atoms are traversed.
    def count_carbons(start_atom, exclude_idxs):
        count = 0
        seen = set()
        queue = [start_atom]
        while queue:
            atom = queue.pop(0)
            if atom.GetIdx() in seen:
                continue
            seen.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:
                count += 1
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in exclude_idxs or nbr.GetIdx() in seen:
                    continue
                # Only traverse through heavy atoms
                if nbr.GetAtomicNum() > 1:
                    queue.append(nbr)
        return count

    # Gather candidate backbone atoms: sp3 carbons not in rings.
    candidate_carbons = [atom for atom in mol.GetAtoms() 
                         if atom.GetAtomicNum() == 6 and 
                            atom.GetHybridization() == rdchem.HybridizationType.SP3 and 
                            not atom.IsInRing()]
    
    # Now try to find a set of three candidate carbons with the correct connectivity:
    # We require a linear chain: two end carbons (c1 and c3) attached to a central carbon (c2)
    # and c1 and c3 not bonded directly.
    for mid in candidate_carbons:
        # Look for neighbors of the mid atom among candidates.
        neighbors = [nbr for nbr in mid.GetNeighbors() 
                     if nbr.GetIdx() in [a.GetIdx() for a in candidate_carbons]]
        if len(neighbors) < 2:
            continue
        # Try unordered pairs of neighbor carbons
        n_count = len(neighbors)
        for i in range(n_count):
            for j in range(i+1, n_count):
                c1 = neighbors[i]
                c3 = neighbors[j]
                # For a true linear chain, the two end carbons should not be directly bonded.
                if mol.GetBondBetweenAtoms(c1.GetIdx(), c3.GetIdx()) is not None:
                    continue
                # We now have a candidate backbone: [c1, mid, c3]
                backbone = [c1, mid, c3]
                backbone_idxs = {atom.GetIdx() for atom in backbone}
                # For each backbone carbon, require exactly one oxygen neighbor not in backbone.
                oxygens_per_carbon = {}
                valid_backbone = True
                for carbon in backbone:
                    oxy_neighbors = [nbr for nbr in carbon.GetNeighbors()
                                     if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_idxs]
                    if len(oxy_neighbors) != 1:
                        valid_backbone = False
                        break
                    oxygens_per_carbon[carbon.GetIdx()] = oxy_neighbors[0]
                if not valid_backbone:
                    continue
                
                # Now for each oxygen attached to the backbone, decide if it is "free" (has H)
                # or "substituted" (no hydrogen but leads to a branch of at least 3 carbons).
                substituted_count = 0
                for carbon in backbone:
                    oxy = oxygens_per_carbon[carbon.GetIdx()]
                    # Use total implicit+explicit H count; if > 0 we assume free –OH.
                    if oxy.GetTotalNumHs() > 0:
                        # free hydroxyl; do nothing
                        continue
                    else:
                        # No hydrogen attached. Check the branch beyond the oxygen.
                        branch_valid = False
                        # Examine oxygen's neighbors other than the backbone carbon.
                        for nbr in oxy.GetNeighbors():
                            if nbr.GetIdx() == carbon.GetIdx():
                                continue
                            # Only consider carbon neighbors
                            if nbr.GetAtomicNum() == 6:
                                # Count the carbons in the branch originating from this neighbor.
                                branch_size = count_carbons(nbr, backbone_idxs.union({oxy.GetIdx()}))
                                if branch_size >= 3:
                                    branch_valid = True
                                    break
                        if branch_valid:
                            substituted_count += 1
                # For a diradylglycerol we need exactly two substituted positions and one free hydroxyl.
                if substituted_count == 2:
                    return True, "Found a glycerol backbone with exactly 2 substituted positions"
    
    return False, "No glycerol backbone with exactly 2 valid substituent branches found"


# Example usage:
if __name__ == "__main__":
    # Test with one example (a DG structure)
    test_smiles = "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"
    result, reason = is_diradylglycerol(test_smiles)
    print("Is Diradylglycerol:", result)
    print("Reason:", reason)