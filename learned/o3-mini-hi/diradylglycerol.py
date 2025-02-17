"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition:
  Any lipid that is glycerol bearing two substituent groups – either acyl, alkyl, or alk-1-enyl –
  at any two of the three possible positions.
  
This improved heuristic:
  1. Parses the SMILES and rejects invalid input.
  2. Immediately excludes molecules containing phosphorus (commonly phospholipids).
  3. Searches for a non‐ring, linear three–carbon “backbone” (candidate glycerol) where each carbon is sp3.
  4. For each backbone carbon, examines non‐backbone oxygen neighbors.
     If an oxygen has an attached hydrogen (using GetTotalNumHs) we interpret it as an unmodified –OH.
     Otherwise (as in an ester or ether), we see if the oxygen leads into a branch containing ≥3 carbons.
  5. If exactly two backbone carbons have a valid substituent branch, we return a positive classification.
  
Note: This is still a heuristic approach and may both reject some valid diradylglycerols and pick up some false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone (a linear chain of three sp3 carbons,
    none in a ring) in which exactly two positions are substituted with an acyl/alkyl/alk-1-enyl chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diradylglycerol, otherwise False.
        str: Explanation for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules that contain phosphorus; these are more likely phospholipids.
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid, not a diradylglycerol"

    # Helper function: count the number of carbon atoms in a connected branch from a given start_atom,
    # while excluding atoms whose indices are in the provided set.
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
                # Traverse only heavy atoms.
                if nbr.GetAtomicNum() > 1:
                    queue.append(nbr)
        return count

    # Find candidate backbone atoms: sp3 carbons that are not in rings.
    candidate_carbons = [atom for atom in mol.GetAtoms() 
                         if atom.GetAtomicNum() == 6 and atom.GetHybridization() == rdchem.HybridizationType.SP3 and not atom.IsInRing()]
    
    # Look for a set of three carbons arranged in a linear chain.
    for mid in candidate_carbons:
        # Get sp3 carbon neighbors (also not in ring) for a candidate middle.
        nbr_carbons = [nbr for nbr in mid.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == rdchem.HybridizationType.SP3 and not nbr.IsInRing()]
        if len(nbr_carbons) < 2:
            continue
        # Try each unordered pair as the two ends.
        for i in range(len(nbr_carbons)):
            for j in range(i+1, len(nbr_carbons)):
                c1 = nbr_carbons[i]
                c3 = nbr_carbons[j]
                # For a linear chain, the two end carbons should not be directly bonded.
                if mol.GetBondBetweenAtoms(c1.GetIdx(), c3.GetIdx()) is not None:
                    continue
                backbone = [c1, mid, c3]
                # Double-check that none of the backbone atoms is in a ring.
                if any(atom.IsInRing() for atom in backbone):
                    continue

                # Now assess the substituents: originally each glycerol carbon bears an –OH.
                # For each backbone carbon, get oxygen neighbors that are not part of the backbone.
                subst_count = 0
                valid_backbone = True
                backbone_idxs = {atom.GetIdx() for atom in backbone}
                for carbon in backbone:
                    oxygens = [nbr for nbr in carbon.GetNeighbors() 
                               if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_idxs]
                    if not oxygens:
                        # If a backbone carbon lacks an oxygen neighbor, it is unlikely to be part of a glycerol structure.
                        valid_backbone = False
                        break
                    # Determine if this carbon is substituted.
                    # A free hydroxyl is expected to have at least one hydrogen on the oxygen.
                    # Otherwise, if the oxygen (with no H attached) leads to a branch with ≥3 carbons,
                    # we consider that position substituted.
                    is_substituted = False
                    for oxy in oxygens:
                        # Check if the oxygen has any hydrogen; if so, it is likely still a free –OH.
                        # (Hydrogens may not be explicit so use GetTotalNumHs.)
                        if oxy.GetTotalNumHs() > 0:
                            continue
                        # Otherwise, check its other neighbors (besides our backbone carbon).
                        for nbr in oxy.GetNeighbors():
                            if nbr.GetIdx() == carbon.GetIdx():
                                continue
                            if nbr.GetAtomicNum() == 6:
                                # Exclude the backbone and this oxygen from the branch search.
                                branch_size = count_carbons(nbr, backbone_idxs.union({oxy.GetIdx()}))
                                if branch_size >= 3:
                                    is_substituted = True
                                    break
                        if is_substituted:
                            break
                    if is_substituted:
                        subst_count += 1
                # For a diradylglycerol, exactly 2 of the 3 backbone carbons must be substituted.
                if valid_backbone and subst_count == 2:
                    return True, "Found a glycerol backbone with exactly 2 substituted positions"
    return False, "No glycerol backbone with exactly 2 valid substituent branches found"

# Example usage:
if __name__ == "__main__":
    # Here is one candidate SMILES for a DG (diradylglycerol) structure.
    test_smiles = "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"
    result, reason = is_diradylglycerol(test_smiles)
    print("Is Diradylglycerol:", result)
    print("Reason:", reason)