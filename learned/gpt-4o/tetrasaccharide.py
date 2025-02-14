"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern to match the furanose and pyranose rings 
    # with variable OH distributions and linkages
    # This pattern is an oversimplified representation and might need refinement
    pyranose_pattern = Chem.MolFromSmarts("[OX2H][C@H]1[C@H]([OX2H])[C@@H]([OX2H])[C@H]([OX2H])[C@@H]1[OX2H]")
    furanose_pattern = Chem.MolFromSmarts("[OX2H][C@H]1[C@H]([OX2H])[C@@H]([OX2H])[C@H]([OX2H])1")
    
    # Count pyranose and furanose rings
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    
    # Each pyranose/furanose ring can be considered a monosaccharide unit
    num_monosaccharide_units = len(pyranose_matches) + len(furanose_matches)
    
    if num_monosaccharide_units == 4:
        return True, "Contains four linked monosaccharide units"
    else:
        return False, f"Found {num_monosaccharide_units} monosaccharide units, need exactly 4"

# Example usage:
# print(is_tetrasaccharide("example_smiles_string"))  # Replace with an actual SMILES for testing