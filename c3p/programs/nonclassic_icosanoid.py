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
    if carbon_count != 20:
        return False, f"Contains {carbon_count} carbons; requires exactly 20 carbons for backbone"
    
    # Check for key functional groups: epoxide, hydroxyl, and carboxylic acid
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")  # Epoxide ring
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # Carboxylate group for stability in aqueous solutions and/or protonated
    
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    carboxyl_count = len(mol.GetSubstructMatches(carboxyl_pattern))

    # Ensure it contains essential functional groups
    if not has_epoxide or hydroxyl_count < 1 or carboxyl_count < 1:
        return False, "Must contain at least one epoxide group, hydroxyl, and carboxylic acid group"
    
    # Avoid classic icosanoid motifs
    prostanoid_pattern = Chem.MolFromSmarts("C1=CCCCC1C(=O)O")  # Simplified prostanoids pattern
    leukotriene_pattern = Chem.MolFromSmarts("CC(C)(C)C")  # Leaf to simplify searching branches

    if mol.HasSubstructMatch(prostanoid_pattern) or mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Substructure resembles classic icosanoids (prostanoid or leukotriene patterns detected)"
    
    # If all criteria are met, classify as nonclassic icosanoid
    return True, "Classified as nonclassic icosanoid with adequate C20 backbone and functionalization"

# Note: This classification focuses on identifying key elements typically found in nonclassic icosanoids 
# while excluding features common to classic icosanoids.