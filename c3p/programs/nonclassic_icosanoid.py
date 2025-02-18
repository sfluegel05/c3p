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
    
    # Step 1: Check for the presence of C20 backbone
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 20:
        return False, f"Contains {carbon_count} carbons; requires exactly 20 carbons"

    # Step 2: Check for presence of functional groups: epoxide, hydroxyl, and carboxylic acid
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)

    # Check if it has at least one of epoxide or hydroxyl group, and a carboxylic acid group
    if not (has_epoxide or has_hydroxyl) or not has_carboxylic_acid:
        return False, "Lacks necessary epoxide or hydroxyl groups, or carboxylic acid group"

    # Step 3: Ensure it's not a classical icosanoid (a more complex task)
    # For simplicity, we may consider this as avoiding common prostanoid and leukotriene motifs
    prostanoid_pattern = Chem.MolFromSmarts("COC(=O)CC")
    leukotriene_pattern = Chem.MolFromSmarts("CCCCC(CC)CC")
    
    if mol.HasSubstructMatch(prostanoid_pattern) or mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Substructure resembles classic icosanoids (prostanoid or leukotriene patterns detected)"
    
    # If all criteria are met, classify as nonclassic icosanoid
    return True, "Classified as nonclassic icosanoid with adequate C20 backbone and functionalization"

# Note: This classification is based on identifying specific structural elements that typically
# feature in nonclassic icosanoids and might not capture all possible derivatives accurately.