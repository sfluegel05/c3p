"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid typically includes a hydroxy-5beta-cholanic acid structure, 
    with specific stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Specific steroid backbone pattern for bile acids with 5beta configuration
    steroid_5beta_pattern = Chem.MolFromSmarts("C1[C@@H]2[C@H](C[C@@H]3[C@H](C2)CC[C@@H]4[C@]3(CCCC4)C)C1")
    if not mol.HasSubstructMatch(steroid_5beta_pattern):
        return False, "Does not match the 5beta steroid backbone of bile acids"
    
    # Check for carboxylic acid group – essential for bile acids
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for multiple hydroxyl groups on the steroid structure
    hydroxyl_pattern = Chem.MolFromSmarts("O[C@H]")  # Commonly present in bile acids
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:  # Bile acids often have multiple hydroxylations
        return False, f"Insufficient hydroxyl groups, found {len(hydroxyl_matches)}"
    
    # Optionally, count chiral centers if needed
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if len(chiral_centers) < 5:
        return False, f"Insufficient chiral centers for bile acid configuration, found {len(chiral_centers)}"
    
    return True, "Matches bile acid structural criteria"