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

    # Define generalized patterns for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H]O[C@H](O)[C@@H](O)[C@H](O)C1")  # 6-membered pyranose
    furanose_pattern = Chem.MolFromSmarts("C1[C@H]O[C@H](O)[C@@H](O)C1")  # 5-membered furanose

    # Define pattern for any glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@@H]")

    # Get all matches for pyranose and furanose patterns
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Calculate total number of rings
    num_sugar_rings = len(pyranose_matches) + len(furanose_matches)

    # Check presence of two sugar rings
    if num_sugar_rings != 2:
        return False, f"Expected 2 sugar rings, found {num_sugar_rings}"

    # Look for glycosidic bond(s)
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    # Ensure there is at least one glycosidic bond
    if not glycosidic_matches:
        return False, "No glycosidic bond found"

    # Double-check if these patterns are correctly connected for classic disaccharide
    # Check if glycosidic oxygen connects two rings
    for gm in glycosidic_matches:
        o_idx = gm[0]
        o_neighbors = mol.GetAtomWithIdx(o_idx).GetNeighbors()

        # Validate that the oxygen connects to 2 sugars
        sugar_count = 0
        for atom in o_neighbors:
            # Check if any of the neighbors of Oxygen are part of sugar rings.
            if atom.GetIdx() in {idx for match in pyranose_matches + furanose_matches for idx in match}:
                sugar_count += 1
        
        if sugar_count == 2:
            return True, "Contains two sugar rings properly connected by a glycosidic bond"

    return False, "Sugar rings are not properly connected by glycosidic bonds"