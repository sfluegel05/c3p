"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid typically includes a hydroxy-5beta-cholanic acid structure, 
    with specific stereochemistry, and can feature amide linkages to glycine or taurine.
    
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
    
    # Relaxed steroid backbone pattern for hydroxy-5beta-cholanic acids
    steroid_5beta_pattern = Chem.MolFromSmarts("C1[C@H]2CC[C@@H]3[C@@]([C@H](CCC4C(CC=1)[C@@]34C)O)C")
    if not mol.HasSubstructMatch(steroid_5beta_pattern):
        return False, "Does not match the hydroxy-5beta-cholanic acid backbone"
    
    # Check for carboxylic acid or carboxamide group
    # Carboxylic acid pattern
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylamide_pattern = Chem.MolFromSmarts("C(=O)N")
    
    if not (mol.HasSubstructMatch(carboxyl_pattern) or mol.HasSubstructMatch(carboxylamide_pattern)):
        return False, "No carboxylic acid or carboxamide group found"
    
    # Ensure the presence of hydroxyl groups, which are common in bile acids
    hydroxyl_pattern = Chem.MolFromSmarts("O[C@H]")  
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxyl_matches)}"
    
    # Validate stereochemistry matches the typical bile acid configuration
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if len(chiral_centers) < 5:
        return False, f"Insufficient chiral centers for bile acid configuration, found {len(chiral_centers)}"
    
    return True, "Matches bile acid structural criteria"