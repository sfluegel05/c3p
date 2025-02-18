"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition:
  Any lipid that is glycerol bearing two substituent groups - either acyl, alkyl, or alk-1-enyl -
  at any two of the three possible positions.
  
Heuristic in this revision:
  1. Parse the SMILES and require a valid molecule.
  2. Exclude molecules with phosphorus (as these usually belong to phospholipids).
  3. Search for an acyclic (non-ring) linear chain of 3 sp3 carbons (candidate glycerol backbone).
     We require that the two end carbons are not directly connected.
  4. For each backbone carbon, collect oxygen neighbors (ignoring backbone atoms).
  5. For each oxygen neighbor (skipping those connected to phosphorus), check if its other neighbor(s)
     (i.e. the branch attached) includes at least 3 carbon atoms in a connected substructure.
  6. If exactly two of the three backbone carbons are substituted with such a branch, then we classify
     the molecule as a diradylglycerol.
     
Note: This is a heuristic approach – it may reject some valid structures or report false positives,
but it is tuned to improve over the previous version.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone (a linear chain of three sp3 carbons,
    none in a ring) with exactly two substituent groups (acyl, alkyl, or alk-1-enyl) attached
    on two of the three positions.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diradylglycerol, otherwise False.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules containing phosphorus; these are typically phospholipids.
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus; likely a phospholipid, not a diradylglycerol"
    
    # Helper function: count the number of carbon atoms reachable from a starting atom,
    # excluding atoms in the given set (e.g. backbone or the oxygen itself).
    def count_carbons(start_atom, exclude_idxs):
        count = 0
        seen = set()
        queue = [start_atom]
        while queue:
            atom = queue.pop(0)
            if atom.GetIdx() in seen:
                continue
            seen.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:  # carbon
                count += 1
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in exclude_idxs or nbr.GetIdx() in seen:
                    continue
                # Only traverse heavy atoms
                if nbr.GetAtomicNum() > 1:
                    queue.append(nbr)
        return count

    # Identify candidate backbone atoms: sp3 carbons that are not in rings.
    candidate_carbons = [atom for atom in mol.GetAtoms() 
                         if atom.GetAtomicNum() == 6 
                         and atom.GetHybridization() == rdchem.HybridizationType.SP3
                         and not atom.IsInRing()]
    
    # Loop over candidate "middle" atoms to build a linear chain backbone
    for mid in candidate_carbons:
        # Get sp3 carbon neighbors (acyclic, not in ring) to serve as possible ends.
        nbr_carbons = [nbr for nbr in mid.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 
                       and nbr.GetHybridization() == rdchem.HybridizationType.SP3 
                       and not nbr.IsInRing()]
        if len(nbr_carbons) < 2:
            continue
        # Try every unordered pair from these neighbors
        for i in range(len(nbr_carbons)):
            for j in range(i+1, len(nbr_carbons)):
                c1 = nbr_carbons[i]
                c3 = nbr_carbons[j]
                # For a linear chain, ensure that the two ends are not directly bonded.
                if mol.GetBondBetweenAtoms(c1.GetIdx(), c3.GetIdx()) is not None:
                    continue
                backbone = [c1, mid, c3]
                # Confirm none of the backbone carbons is in a ring.
                if any(atom.IsInRing() for atom in backbone):
                    continue

                # Now, for each backbone carbon, assess substitution.
                # We expect that in a genuine glycerol backbone, each carbon was originally –OH substituted.
                # Here, a substituent is considered valid if the oxygen attached to the carbon (and not part of the backbone)
                # is connected to a branch that has at least 3 carbon atoms (heuristic for an acyl/alkyl chain).
                substituted_count = 0
                valid_backbone = True
                backbone_idxs = {a.GetIdx() for a in backbone}
                for carbon in backbone:
                    # Get oxygen neighbors not belonging to the backbone.
                    oxy_neighbors = [nbr for nbr in carbon.GetNeighbors() 
                                     if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_idxs]
                    if not oxy_neighbors:
                        # A backbone carbon without any oxygen is unlikely to be a glycerol-derived carbon.
                        valid_backbone = False
                        break
                    # For this carbon, see if at least one oxygen carries a sufficiently long substituent.
                    found_valid_subs = False
                    for oxy in oxy_neighbors:
                        # Skip an oxygen that is connected to phosphorus (e.g. in phospholipids).
                        if any(nbr.GetAtomicNum() == 15 for nbr in oxy.GetNeighbors()):
                            continue
                        # Examine each neighbor of the oxygen, excluding the backbone carbon.
                        for nbr in oxy.GetNeighbors():
                            if nbr.GetIdx() == carbon.GetIdx():
                                continue
                            # Check if neighbor is carbon and its branch (excluding oxy and backbone)
                            if nbr.GetAtomicNum() == 6:
                                exclude_set = backbone_idxs.union({oxy.GetIdx()})
                                branch_size = count_carbons(nbr, exclude_set)
                                if branch_size >= 3:
                                    found_valid_subs = True
                                    break
                        if found_valid_subs:
                            break
                    if found_valid_subs:
                        substituted_count += 1
                # For a diradylglycerol we require exactly two substituted positions.
                if valid_backbone and substituted_count == 2:
                    return True, "Found a glycerol backbone with exactly 2 substituted positions"
    return False, "No glycerol backbone found with exactly 2 substituent groups"

# Example usage: (this part could be run to test the classification)
if __name__ == "__main__":
    # A test SMILES string for a DG (diradylglycerol) candidate:
    test_smiles = "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"
    result, reason = is_diradylglycerol(test_smiles)
    print("Is Diradylglycerol:", result)
    print("Reason:", reason)