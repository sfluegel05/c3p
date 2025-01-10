"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a tetrasaccharide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for monosaccharide units - pyranose (6-membered) and furanose (5-membered) rings
    pyranose_pattern = Chem.MolFromSmarts("[C@H1]([O])[C@H1][C2@H][C@H2][C@H2][C@H2]([O])O2") # Considering typical pyranose rings
    furanose_pattern = Chem.MolFromSmarts("[C@H]1[O][C@@H][C][C@H1]") # Considering typical furanose rings
    
    # Count matches for pyranose and furanose rings
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Each ring generally signifies a monosaccharide unit in saccharide chemistry
    num_monosaccharide_units = len(pyranose_matches) + len(furanose_matches)
    
    if num_monosaccharide_units == 4:
        return True, "Contains four linked monosaccharide units"
    else:
        return False, f"Found {num_monosaccharide_units} monosaccharide units, need exactly 4"

# Example usage - Replace with actual SMILES for testing
# print(is_tetrasaccharide("some_smiles_string"))