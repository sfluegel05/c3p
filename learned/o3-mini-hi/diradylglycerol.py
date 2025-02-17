"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition:
  Any lipid having a glycerol backbone (a linear chain of three sp3 carbons, none in a ring)
  bearing exactly two substituent groups (acyl, alkyl, or alk-1-enyl) on two of the 
  three possible positions.
  
Heuristic in this revision:
  1. Parse the SMILES and require a valid molecule.
  2. Search for a noncyclic (acyclic) linear chain of three sp³ carbons 
     (candidate backbone: A-B-C where B is connected to A and C and A and C are not bonded).
  3. For each backbone carbon, obtain oxygen neighbors (excluding those that are part of the backbone).
  4. For each oxygen neighbor, check if it is “substituted.” We do this by looking at its 
     non-backbone neighbors: if one of them is a carbon and the connected carbon substructure 
     (found by a small breadth-first search) has at least 3 carbon atoms, we mark that backbone position 
     as substituted.
  5. If exactly two of the three backbone carbons are substituted, then classify as diradylglycerol.
  
This heuristic is tuned to reduce false positives from cyclic or highly decorated molecules
and to capture genuine diradylglycerols.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone (a linear chain of three sp3 carbons,
    none in a ring) with exactly two substituent groups (acyl, alkyl, or alk-1-enyl)
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
    
    # Helper: count the number of carbon atoms in the substituent branch reachable
    # from a given neighbor (excluding atoms in the 'exclude_idxs' set).
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
                # Only traverse heavy atoms (ignore hydrogens which are implicit usually)
                if nbr.GetAtomicNum() > 1:
                    queue.append(nbr)
        return count
    
    # Get all sp3 carbons that are not in rings (we want acyclic backbone)
    sp3_carbons = [atom for atom in mol.GetAtoms() 
                   if atom.GetAtomicNum() == 6 
                   and atom.GetHybridization() == rdchem.HybridizationType.SP3
                   and not atom.IsInRing()]
    
    # For each sp3 carbon, try to use it as the middle of a linear chain (backbone candidate)
    for mid in sp3_carbons:
        # Get sp3 carbon neighbors that are acyclic and not in a ring
        carbon_neighbors = [nbr for nbr in mid.GetNeighbors() 
                            if nbr.GetAtomicNum() == 6 
                            and nbr.GetHybridization() == rdchem.HybridizationType.SP3 
                            and not nbr.IsInRing()]
        if len(carbon_neighbors) < 2:
            continue
        # Try every unordered pair from these neighbors
        for i in range(len(carbon_neighbors)):
            for j in range(i+1, len(carbon_neighbors)):
                c1 = carbon_neighbors[i]
                c3 = carbon_neighbors[j]
                # Ensure that c1 and c3 are not directly bonded (to keep it linear, not cyclic)
                if mol.GetBondBetweenAtoms(c1.GetIdx(), c3.GetIdx()) is not None:
                    continue
                backbone = [c1, mid, c3]
                # For good measure, re-check that these three atoms are distinct and acyclic.
                if any(atom.IsInRing() for atom in backbone):
                    continue
                
                # Now for each backbone carbon, look for oxygen neighbors that are not in
                # the backbone. We expect a true glycerol carbon to have at least one oxygen (free or substituted).
                substituted_count = 0
                valid_backbone = True
                for carbon in backbone:
                    # Find oxygen neighbors not in backbone.
                    backbone_idxs = {atom.GetIdx() for atom in backbone}
                    oxy_neighbors = [nbr for nbr in carbon.GetNeighbors() 
                                     if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_idxs]
                    if not oxy_neighbors:
                        # If a backbone carbon lacks any oxygen neighbor, this candidate is unlikely glycerol.
                        valid_backbone = False
                        break
                    # Check if at least one oxygen neighbor carries a substituent.
                    found_substituent = False
                    for oxy in oxy_neighbors:
                        # For each oxygen bound to the backbone carbon, examine its other neighbors.
                        # Exclude the backbone carbon.
                        for nbr in oxy.GetNeighbors():
                            if nbr.GetIdx() == carbon.GetIdx():
                                continue
                            # We consider a substituent valid if it is a carbon and is part
                            # of a branch that has at least 3 carbons (heuristic for acyl/alkyl chain).
                            if nbr.GetAtomicNum() == 6:
                                # Do a BFS starting at this neighbor; exclude the oxygen and backbone atoms.
                                exclude_set = backbone_idxs.union({oxy.GetIdx()})
                                chain_size = count_carbons(nbr, exclude_set)
                                if chain_size >= 3:
                                    found_substituent = True
                                    break
                        if found_substituent:
                            break
                    if found_substituent:
                        substituted_count += 1
                # We require that exactly two of the three backbone carbons are substituted
                if valid_backbone and substituted_count == 2:
                    return True, "Found a glycerol backbone with exactly 2 substituted positions"
    return False, "No glycerol backbone found with exactly 2 substituent groups"

# Example usage:
if __name__ == "__main__":
    # Test with an example diradylglycerol SMILES
    test_smiles = "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"
    result, reason = is_diradylglycerol(test_smiles)
    print("Is Diradylglycerol:", result)
    print("Reason:", reason)