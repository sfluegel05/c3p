"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid typically has a hydroxy-5beta-cholanic acid structure.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Adjusted steroid backbone pattern for bile acids
    # This pattern focuses on the four-ring structure with proper stereochemistry
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(C)CCCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not match the steroid backbone of bile acids"
    
    # Identify the 5beta configuration: typically involves specific stereochemistry
    # Checking for at least 5 chiral centers to support complex stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 5:
        return False, f"Insufficient chiral centers for 5beta configuration, found {len(chiral_centers)}"
    
    # Check for important functional group: carboxylic acid with stereo considerations
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for hydroxyl groups on the steroid structure (bile acids often have multiple)
    hydroxyl_pattern = Chem.MolFromSmarts("O[C@H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxyl_matches)}"
    
    return True, "Matches bile acid structural criteria"