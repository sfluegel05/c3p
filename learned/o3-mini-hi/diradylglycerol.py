"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition:
  Any lipid having a glycerol backbone (a linear chain of 3 sp3 carbons) 
  bearing two substituent groups (acyl, alkyl, or alk-1-enyl) at any two of the 
  three possible positions.
  
Heuristic approach in this revision:
  1. Parse the SMILES string.
  2. Identify candidate glycerol backbones by finding a linear (noncyclic) chain 
     of three sp³ carbons. (We require that the two terminal carbons are not directly
     connected.)
  3. For each candidate backbone carbon (excluding the bond joining the three), 
     search its oxygen neighbors (excluding those that belong to the backbone).
  4. For each oxygen neighbor we consider it “substituted” if it is attached to at least 
     two heavy (atomic num > 1) atoms – one being the backbone carbon and at least one additional.
  5. If exactly two of the three backbone carbons have at least one such substituted oxygen,
     then we classify the molecule as a diradylglycerol.
  
If no candidate backbone meets these criteria, the function returns False with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone (a linear chain of three sp3 carbons)
    bearing exactly two substituent groups (acyl, alkyl, or alk-1-enyl) on two of its three positions.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diradylglycerol, otherwise False.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We restrict ourselves to sp3 carbons
    sp3_carbons = [atom for atom in mol.GetAtoms() 
                   if atom.GetAtomicNum() == 6 
                   and atom.GetHybridization() == rdchem.HybridizationType.SP3]
    
    # Look for a 3-atom linear chain: candidate backbone = (C1, C2, C3)
    # where C2 (middle) is connected to both C1 and C3, and C1 and C3 are not directly bonded.
    for mid in sp3_carbons:
        # Get sp3 carbon neighbors of mid (candidates for terminal positions)
        carbon_neighbors = [nbr for nbr in mid.GetNeighbors() 
                            if nbr.GetAtomicNum() == 6 
                            and nbr.GetHybridization() == rdchem.HybridizationType.SP3]
        # We need at least two; try every unordered pair.
        if len(carbon_neighbors) < 2:
            continue
        for i in range(len(carbon_neighbors)):
            for j in range(i+1, len(carbon_neighbors)):
                c1 = carbon_neighbors[i]
                c3 = carbon_neighbors[j]
                # Ensure c1 and c3 are not directly bonded (to avoid cycles)
                if mol.GetBondBetweenAtoms(c1.GetIdx(), c3.GetIdx()) is not None:
                    continue
                backbone = [c1, mid, c3]
                # For each backbone carbon, search for oxygen neighbors that are not in the backbone.
                # We mark the carbon as having a substituent if at least one such oxygen is “substituted.”
                substituted_positions = 0
                valid_backbone = True
                for carbon in backbone:
                    # Get oxygen neighbors not in the backbone.
                    oxy_neighbors = [nbr for nbr in carbon.GetNeighbors() 
                                     if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [a.GetIdx() for a in backbone]]
                    # For a glycerol, we expect each backbone carbon to have at least one oxygen.
                    if not oxy_neighbors:
                        valid_backbone = False
                        break
                    # Check if any oxygen neighbor qualifies as substituted.
                    found_substituted = False
                    for oxy in oxy_neighbors:
                        # Count heavy neighbors (atomic number > 1). 
                        # A free hydroxyl oxygen normally bonds only to the backbone (and an implicit hydrogen).
                        heavy_count = sum(1 for n in oxy.GetNeighbors() if n.GetAtomicNum() > 1)
                        if heavy_count > 1:
                            found_substituted = True
                            break
                    if found_substituted:
                        substituted_positions += 1
                # We expect exactly two of the three backbone carbons to have a substituted oxygen.
                if valid_backbone and substituted_positions == 2:
                    return True, "Found a glycerol backbone with exactly 2 substituted positions"
    return False, "No glycerol backbone found with exactly 2 substituent groups"


# Example usage:
if __name__ == "__main__":
    # We test with one provided SMILES (a diradylglycerol example)
    test_smiles = "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"
    result, reason = is_diradylglycerol(test_smiles)
    print("Is Diradylglycerol:", result)
    print("Reason:", reason)