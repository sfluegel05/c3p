"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Identify ketal carbons: carbon bonded to two oxygens (sp3 hybridized) and two carbons, no hydrogens
    ketal_pattern = Chem.MolFromSmarts("[CX4](-[OX2])(-[OX2])([CX4])([CX4])") # Modified to explicitly request no hydrogens.
    ketal_matches = mol.GetSubstructMatches(ketal_pattern)
    if not ketal_matches:
        return False, "No ketal carbon found."

    # Iterate through each ketal carbon match
    for match in ketal_matches:
      ketal_carbon_idx = match[0]
      ketal_carbon = mol.GetAtomWithIdx(ketal_carbon_idx)

      # Get oxygen atoms bonded to ketal carbon
      oxygen_indices = [neighbor.GetIdx() for neighbor in ketal_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 8]
      
      if len(oxygen_indices) != 2:
        continue # Should not happen because of the pattern but just in case.
      
      oxy1_idx = oxygen_indices[0]
      oxy2_idx = oxygen_indices[1]
      oxy1 = mol.GetAtomWithIdx(oxy1_idx)
      oxy2 = mol.GetAtomWithIdx(oxy2_idx)

      # 2. Verify that both oxygens are part of rings. 
      if not oxy1.IsInRing() or not oxy2.IsInRing():
            continue

      # 3. Verify if the two rings share only the ketal carbon atom
      ring_info = mol.GetRingInfo()

      rings_oxy1 = [ring for ring in ring_info.AtomRings() if oxy1_idx in ring]
      rings_oxy2 = [ring for ring in ring_info.AtomRings() if oxy2_idx in ring]
      
      if not rings_oxy1 or not rings_oxy2:
          continue
      
      # Select the first ring that each oxygen is part of
      ring1 = rings_oxy1[0]
      ring2 = rings_oxy2[0]

      # If rings are the same, not a spiroketal
      if set(ring1) == set(ring2):
        continue
    
      # Check that intersection of the two rings are only the spiroketal carbon
      common_atoms = set(ring1).intersection(set(ring2))
      if len(common_atoms) == 1 and ketal_carbon_idx in common_atoms:
        return True, "Spiroketal structure identified."
    
    return False, "No spiroketal structure detected."