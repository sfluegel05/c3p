"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharide units linked by a glycosidic bond.

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

    # Define potential patterns for sugar rings (pyranose and furanose) and glycosidic bond
    pyranose_pattern = Chem.MolFromSmarts("C1OC([H,O])([H,O])C([H,O])([H,O])C([H,O])([H,O])O1")
    furanose_pattern = Chem.MolFromSmarts("C1OC([H,O])([H,O])C([H,O])([H,O])O1")
    glycosidic_bond_pattern = Chem.MolFromSmarts("[O]C(O)O")
    
    # Determine the number of sugar rings
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    num_sugar_rings = len(pyranose_matches) + len(furanose_matches)

    if num_sugar_rings != 2:
        return False, f"Expected 2 sugar rings, found {num_sugar_rings}"

    # Check for glycosidic bonds
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    if not glycosidic_matches:
        return False, "No glycosidic bond found"

    return True, "Contains two sugar rings linked by a glycosidic bond"