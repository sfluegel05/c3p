"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic Icosanoid
"""
from rdkit import Chem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of a central C20 backbone
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Contains {carbon_count} carbons; requires at least 20 carbons for backbone"
    
    # Check for presence of functional groups: epoxide, multiple hydroxyls, and carboxylic acid
    epoxide_pattern = Chem.MolFromSmarts("[CX3]1[OX2][CX3]1")  # Improved pattern for epoxides
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Matches hydroxyl groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")  # Matches carboxylic acid group
    
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern)) >= 2  # Require two or more hydroxyl groups
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)

    # Check if it has necessary functional groups
    if not has_epoxide or not hydroxyl_matches or not has_carboxylic_acid:
        return False, "Must contain epoxide group, multiple hydroxyls, and a carboxylic acid group"
    
    # Avoid classic icosanoid motifs
    # Updated SMARTS for more specific avoidance in classical icosanoids
    prostanoid_pattern = Chem.MolFromSmarts("[CX3](=O)[O][CX2][CX3](=O)[CX4]")
    leukotriene_pattern = Chem.MolFromSmarts("CC(C)C[CX4]")
    
    if mol.HasSubstructMatch(prostanoid_pattern) or mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Substructure resembles classic icosanoids (prostanoid or leukotriene patterns detected)"
    
    # If all criteria are met, classify as nonclassic icosanoid
    return True, "Classified as nonclassic icosanoid with adequate C20 backbone and functionalization"

# Note: This classification is based on identifying specific structural elements that typically
# feature in nonclassic icosanoids and might not capture all possible derivatives accurately.