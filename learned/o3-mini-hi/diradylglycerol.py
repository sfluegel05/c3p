"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: Diradylglycerol
Definition:
  Any lipid that is a glycerol bearing two substituent groups – either acyl, alkyl,
  or alk-1-enyl – at any two of the three possible positions.
  
Heuristic approach:
  1. Parse the SMILES string.
  2. Loop through candidate sets of three carbons that are connected (a potential glycerol backbone).
  3. For each candidate backbone [C1, C2, C3] (with C2 connected to both C1 and C3)
     check that each carbon has at least one oxygen neighbor.
  4. For each oxygen attached to a backbone carbon, decide if it is “substituted.”
     (Here, an oxygen is considered substituted if it is connected to two heavy atoms;
      a free hydroxyl oxygen normally only connects to the backbone carbon because its hydrogen is implicit.)
  5. If exactly two of the three backbone carbons show a substituted oxygen, we classify the molecule as diradylglycerol.
  
If no appropriate candidate fingerprint is found, the function returns False with a reason.
"""

from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone bearing exactly two substituent groups
    (acyl, alkyl, or alk-1-enyl) at any two of the three possible positions.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diradylglycerol, otherwise False.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic: Try to locate a glycerol backbone.
    # We look for a chain of 3 carbons (candidate backbone: C1 - C2 - C3) such that each carbon has at least one attached oxygen.
    # Then, for each backbone carbon, we check its oxygen neighbors – if the oxygen in question is connected
    # to another heavy atom (atomic number > 1) besides the backbone carbon then we count it as a substituent.
    
    # Get all carbon atoms in the molecule.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Loop over carbons to pick one as the middle carbon (C2) and then try pairs of its carbon neighbors.
    for mid in carbons:
        # Only consider atoms that have at least two carbon neighbors.
        c_neighbors = [nbr for nbr in mid.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(c_neighbors) < 2:
            continue
        # Try each unordered pair of distinct neighbors as candidate terminal carbons (C1 and C3).
        for i in range(len(c_neighbors)):
            for j in range(i+1, len(c_neighbors)):
                c1 = c_neighbors[i]
                c3 = c_neighbors[j]
                # We now have candidate backbone [c1, mid, c3]
                candidate = [c1, mid, c3]
                valid_backbone = True
                substituent_count = 0
                # For each carbon in the candidate backbone, assess oxygen attachments.
                for c in candidate:
                    # Gather oxygen neighbors that are NOT part of the backbone.
                    oxy_neighbors = [nbr for nbr in c.GetNeighbors() if nbr.GetAtomicNum() == 8]
                    if not oxy_neighbors:
                        valid_backbone = False
                        break  # This candidate does not resemble a glycerol center.
                    # For this backbone carbon, check if at least one oxygen appears substituted.
                    # (Assume that free hydroxyl oxygen has only one heavy (i.e. non-hydrogen) neighbor,
                    #  while an oxygen forming an ester/ether bond will have two heavy neighbors.)
                    found_substituted = False
                    for oxy in oxy_neighbors:
                        heavy_neighbors = [nbr for nbr in oxy.GetNeighbors() if nbr.GetAtomicNum() > 1]
                        # If oxygen is connected to more than one heavy atom then it is substituted.
                        # (Note: the implicit hydrogens are not counted since they are not explicit heavy atoms.)
                        if len(heavy_neighbors) > 1:
                            found_substituted = True
                            break
                    if found_substituted:
                        substituent_count += 1
                if not valid_backbone:
                    continue
                # For diradylglycerol we expect exactly two of the three backbone oxygens to be substituted.
                if substituent_count == 2:
                    return True, "Found a glycerol backbone with exactly 2 substituted positions"
    return False, "No glycerol backbone found with exactly 2 substituent groups"

# Example usage:
if __name__ == "__main__":
    # test with one of the provided SMILES for a DG (diradylglycerol)
    test_smiles = "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"
    result, reason = is_diradylglycerol(test_smiles)
    print("Is Diradylglycerol:", result)
    print("Reason:", reason)