"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
A D-galactoside having beta-configuration at its anomeric centre
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_beta_D_galactoside, reason)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit Hs - this helps with stereochemistry
    mol = Chem.AddHs(mol)
    
    # Beta-D-galactopyranoside core pattern with stereochemistry
    # [OH1] represents the glycosidic oxygen
    # First carbon is the anomeric carbon with beta configuration
    # Fourth carbon has axial OH (characteristic of galactose)
    beta_gal_pattern = Chem.MolFromSmarts(
        "[OD2]-[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H]([CH2]O)O1"
    )
    
    if beta_gal_pattern is None:
        return None, "Invalid SMARTS pattern"
    
    # Look for matches using stereochemistry
    params = Chem.SubstructMatchParameters()
    params.useChirality = True
    matches = mol.GetSubstructMatches(beta_gal_pattern, params)
    
    if not matches:
        return False, "No beta-D-galactoside moiety found"
        
    # Verify glycosidic bond
    for match in matches:
        glycosidic_O = mol.GetAtomWithIdx(match[0])
        
        # Check if oxygen is part of a glycosidic bond (connected to two carbons)
        carbon_neighbors = [n for n in glycosidic_O.GetNeighbors() 
                          if n.GetAtomicNum() == 6]
        
        if len(carbon_neighbors) == 2:
            # Found a valid glycosidic bond
            return True, "Contains beta-D-galactoside moiety with correct stereochemistry"
    
    return False, "Found galactose ring but no glycosidic bond"

def test_smiles():
    """Test function with a known beta-D-galactoside"""
    test_mol = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"  # methyl beta-D-galactoside
    result, reason = is_beta_D_galactoside(test_mol)
    print(f"Test result: {result}, Reason: {reason}")

if __name__ == "__main__":
    test_smiles()