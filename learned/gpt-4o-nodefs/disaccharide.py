"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two sugar rings connected via a proper glycosidic bond.

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

    # Define patterns for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts("C1OC([CH2])[CH2]OC(C1)O")  # 6-membered pyranose pattern
    furanose_pattern = Chem.MolFromSmarts("C1OC([CH2])[CH2]OC(1)O")   # 5-membered furanose pattern

    # Define pattern for glycosidic bond (O connecting two carbon rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C&!R][C&!R]")     

    # Get all pyranose and furanose matches
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Calculate total number of sugar rings
    num_sugar_rings = len(pyranose_matches) + len(furanose_matches)

    # Check for the presence of at least two sugar rings
    if num_sugar_rings < 2:
        return False, f"Expected at least 2 sugar rings, found {num_sugar_rings}"

    # Check for glycosidic bonds
    glycosidic_match_count = len(mol.GetSubstructMatches(glycosidic_bond_pattern))
    if glycosidic_match_count == 0:
        return False, "No glycosidic bond found"

    # Ensure sugar rings are connected via glycosidic bond(s)
    atom_pairs = [(match[0], match[1]) for match in mol.GetSubstructMatches(glycosidic_bond_pattern)]
    for (a1, a2) in atom_pairs:
        for (r1, r2) in pyranose_matches + furanose_matches:
            if a1 in (r1, r2) or a2 in (r1, r2):
                return True, "Contains two sugar rings connected by a glycosidic bond"
            
    return False, "Sugar rings are not properly connected by glycosidic bonds"