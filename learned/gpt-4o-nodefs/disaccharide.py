"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two sugar rings connected via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved patterns for sugar rings and glycosidic bond
    pyranose_pattern = Chem.MolFromSmarts("C1OC([H,O])C([H,O])C([H,O])C([H,O])O1") # 6-membered ring
    furanose_pattern = Chem.MolFromSmarts("C1OC([H,O])C([H,O])C([H,O])O1") # 5-membered ring
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2H]C([OX2H])[!#6]") # Simplified glycosidic bond
    
    # Get all the pyranose and furanose ring matches
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Calculate the total number of sugar rings
    num_sugar_rings = len(pyranose_matches) + len(furanose_matches)
    
    # Check for the presence of at least two sugar rings
    if num_sugar_rings < 2:
        return False, f"Expected at least 2 sugar rings, found {num_sugar_rings}"

    # Check for glycosidic bonds
    glycosidic_match = mol.HasSubstructMatch(glycosidic_bond_pattern)

    if not glycosidic_match:
        return False, "No glycosidic bond found"

    return True, "Contains two sugar rings connected by a glycosidic bond"