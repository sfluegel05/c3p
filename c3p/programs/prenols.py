"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: Prenols – any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH 
in which the carbon skeleton is composed of one or more isoprene units.

This improved function uses a heuristic approach:
  1. It first checks that the SMILES string represents a valid, acyclic molecule.
  2. It then confirms that at least one free hydroxyl group is present.
     (It excludes hydroxyl groups that are part of carboxylic acids.)
  3. It uses two SMARTS definitions so as to detect an isoprene unit regardless of bond direction.
  4. It also checks that there is a minimum number of carbon atoms.
  
If all these conditions are met, the molecule is classified as a prenol.
"""

from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    
    A prenol is defined as an alcohol with a carbon skeleton built exclusively from
    one or more isoprene units. This heuristic checks that:
      - The molecule is valid and acyclic.
      - It contains at least one free hydroxyl (-OH) group (and not an acid –OH).
      - It contains at least one isoprene unit. For this purpose, two SMARTS are applied:
          • isoprene_forward: [CH2]=[C]([CH3])[CH2]
          • isoprene_reverse: [CH2][C]([CH3])=[CH2]
      - It has at least 5 carbon atoms (e.g. isoprenol).
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule is a prenol, False otherwise.
      str: Explanation of the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject cyclic molecules (prenols are acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic, not a prenol"
    
    # Check for the presence of any alcohol group: [OX2H]
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if not alcohol_matches:
        return False, "No hydroxyl (alcohol) group found"
    
    # Exclude hydroxyl groups that are part of carboxylic acids.
    # Carboxylic acid OH groups match the SMARTS C(=O)[OX2H]
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    # If all found hydroxyls are part of a carboxyl group, then no free alcohol is present.
    if len(alcohol_matches) <= len(carboxyl_matches):
        return False, "Hydroxyl group found only in carboxylic acid; no free alcohol present"
    
    # Define two SMARTS patterns to capture an isoprene unit.
    # The isoprene unit here is defined as the fragment:
    #   CH2=C(CH3)CH2   (forward direction) or CH2C(CH3)=CH2 (reverse direction)
    isoprene_forward = Chem.MolFromSmarts("[CH2]=[C]([CH3])[CH2]")
    isoprene_reverse = Chem.MolFromSmarts("[CH2][C]([CH3])=[CH2]")
    
    matches_forward = mol.GetSubstructMatches(isoprene_forward)
    matches_reverse = mol.GetSubstructMatches(isoprene_reverse)
    total_isoprene = len(matches_forward) + len(matches_reverse)
    
    if total_isoprene < 1:
        return False, "No isoprene unit ([CH2]=[C]([CH3])[CH2] or reverse) found"
    
    # Check that the molecule has a minimal number of carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, "Too few carbon atoms to be a prenol"
    
    # All heuristic checks passed
    return True, f"Contains {total_isoprene} isoprene unit(s) and a free hydroxyl group, is acyclic and has {carbon_count} carbons"

# Example usage:
if __name__ == "__main__":
    # Several example SMILES strings can be tested; here is one simple prenol (isoprenol):
    test_smiles = "CC(C)=CCO"  
    result, reason = is_prenols(test_smiles)
    print("Result:", result)
    print("Reason:", reason)