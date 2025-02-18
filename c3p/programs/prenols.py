"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: Prenols – any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH 
in which the carbon skeleton is exclusively composed of one or more isoprene units.

This improved function uses a heuristic approach:
  1. Ensure the SMILES string parses into a valid, acyclic molecule.
  2. Check that at least one free (non–carboxylic) hydroxyl group (-OH) is present.
  3. Exclude molecules that contain phosphorus (such as diphosphates) or other non-terpenoid elements.
  4. Look for evidence in the structure of an isoprene unit via two SMARTS patterns:
       • Pattern1: [CH3]-[C]([CH3])=[CH]-[CH2]
       • Pattern2: [CH2]-[C]([CH3])=[CH]-[CH2]
  5. Also check that there is a minimal number of carbons.
  
If all these conditions are met, the molecule is classified as a prenol.
"""

from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    
    A prenol is defined as an alcohol whose carbon skeleton is 
    composed exclusively of one or more isoprene units.
    
    Heuristic steps:
      - The molecule must be valid and acyclic.
      - It must contain at least one free hydroxyl (alcohol) group 
        (i.e. not all hydroxyl groups are part of a carboxyl group).
      - It must not contain phosphorus atoms (which are seen in diphosphates).
      - It must contain at least one isoprene unit. For our purposes, two SMARTS patterns are used:
          Pattern1 (chain starts with CH3): "[CH3]-[C]([CH3])=[CH]-[CH2]"
          Pattern2 (chain starts with CH2): "[CH2]-[C]([CH3])=[CH]-[CH2]"
      - It must have at least 5 carbon atoms.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule is classified as a prenol, False otherwise.
      str: Explanation of the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject cyclic molecules
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic, not a prenol"
    
    # Exclude molecules that contain phosphorus (e.g., diphosphates)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus, not a prenol"
    
    # Check for the presence of any hydroxyl (-OH) group.
    alcohol_smarts = "[OX2H]"
    alcohol_pattern = Chem.MolFromSmarts(alcohol_smarts)
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if not alcohol_matches:
        return False, "No hydroxyl (alcohol) group found"
    
    # Exclude -OH groups that belong to carboxylic acids.
    carboxyl_smarts = "C(=O)[OX2H]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    # If all hydroxyls are part of carboxylic acid groups, then there is no free alcohol.
    if len(alcohol_matches) <= len(carboxyl_matches):
        return False, "Alcohol group found only in carboxyl functionality; no free hydroxyl present"
    
    # Define two SMARTS patterns to capture an isoprene unit.
    # Pattern1: chain starting with CH3: CH3-C(CH3)=CH-CH2
    isoprene_pat1 = Chem.MolFromSmarts("[CH3]-[C]([CH3])=[CH]-[CH2]")
    # Pattern2: chain starting with CH2: CH2-C(CH3)=CH-CH2
    isoprene_pat2 = Chem.MolFromSmarts("[CH2]-[C]([CH3])=[CH]-[CH2]")
    
    matches_pat1 = mol.GetSubstructMatches(isoprene_pat1)
    matches_pat2 = mol.GetSubstructMatches(isoprene_pat2)
    total_isoprene = len(matches_pat1) + len(matches_pat2)
    
    if total_isoprene < 1:
        return False, "No isoprene unit (either [CH3]-C([CH3])=[CH]-[CH2] or [CH2]-C([CH3])=[CH]-[CH2]) found"
    
    # Count the total number of carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, f"Too few carbon atoms ({carbon_count}) to be a prenol"
    
    # Heuristic: In a prenol made of n isoprene units, roughly 4*n+? carbons are expected.
    # If there are far too many carbons per isoprene unit, it might indicate extra functionality.
    if carbon_count > total_isoprene * 8:
        return False, f"Unexpected high number of carbons ({carbon_count}) relative to detected isoprene unit(s) ({total_isoprene})"
    
    return True, (f"Classified as prenol: contains {total_isoprene} isoprene unit(s), " +
                  f"{len(alcohol_matches) - len(carboxyl_matches)} free hydroxyl(s) and {carbon_count} carbons.")

# Example usage:
if __name__ == "__main__":
    # Test using a simple prenol (isoprenol) example.
    test_smiles = "CC(C)=CCO"  # (2-cis,6-trans)-farnesol analog (mini version)
    result, reason = is_prenols(test_smiles)
    print("Result:", result)
    print("Reason:", reason)