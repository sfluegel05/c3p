"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: Prenols – any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH 
           in which the carbon skeleton is composed of one or more isoprene units.
           
This function uses RDKit to check for a free hydroxyl group and at least one isoprene unit.
The isoprene unit is approximated by the SMARTS pattern "[CH2]=[C]([CH3])[CH2]" which 
captures the characteristic fragment present in many prenol structures.
"""

from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    
    A prenol is an alcohol that incorporates one or more isoprene units.
    This heuristic checks that:
      1. The molecule is valid.
      2. It contains at least one free hydroxyl (-OH) group.
      3. It contains at least one isoprene unit, approximated by the pattern [CH2]=[C]([CH3])[CH2].
      4. It contains at least one C=C double bond.
      5. It has a minimum number of carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
      
    Returns:
        bool: True if the molecule meets the prenol criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of an alcohol group (hydroxyl; [OX2H])
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No free hydroxyl (alcohol) group found"
    
    # Check for the isoprene unit.
    # Here we define an isoprene unit as the fragment [CH2]=[C]([CH3])[CH2]
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C]([CH3])[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if not isoprene_matches:
        return False, "No isoprene unit ([CH2]=[C]([CH3])[CH2]) found"
    
    # Check that at least one non-aromatic C=C double bond is present.
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No carbon–carbon double bond found"
    
    # Optionally, check that the molecule has a minimal number of carbon atoms.
    # Prenols (even the simplest isoprenol) require at least 5 carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, "Too few carbon atoms to be considered a prenol"
    
    return True, f"Contains {len(isoprene_matches)} isoprene unit(s) and a free hydroxyl group"

# Example usage:
if __name__ == "__main__":
    test_smiles = "CC(C)=CCO"  # example prenol SMILES (isoprenol)
    result, reason = is_prenols(test_smiles)
    print("Result:", result)
    print("Reason:", reason)