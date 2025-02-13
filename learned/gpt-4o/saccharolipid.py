"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid contains a carbohydrate moiety linked to a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a pattern for a carbohydrate moiety: a cyclic structure with oxygen atoms (e.g., monosaccharides)
    carbohydrate_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C(O)C1")  # Hexopyranose ring pattern
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No carbohydrate moiety pattern found"

    # Define a pattern for a lipid moiety: hydrocarbon chain with an ester/amide linkage
    lipid_chain_pattern = Chem.MolFromSmarts("C(=O)[O,N][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")  # Fatty acid ester/amide
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No lipid moiety pattern found"

    # Check for any plausible linkage between the carbohydrate and lipid components
    carb_matches = mol.GetSubstructMatches(carbohydrate_pattern)
    lipid_matches = mol.GetSubstructMatches(lipid_chain_pattern)

    # A simplistic connectivity check assuming spatial proximity and shared structure linkages for saccharolipids
    for carb_match in carb_matches:
        for lipid_match in lipid_matches:
            # Simplifies to assume connectedness without explicit bond data, as rdkit lacks direct linkage detection between substructures
            if set(carb_match).intersection(lipid_match):
                return True, "Contains both carbohydrate and lipid components directly linked"

    return False, "Carbohydrate and lipid components found but not linked within the molecule"

# Note: Adjust the SMARTS patterns and connectivity logic as needed since chemistry varies widely among saccharolipids.