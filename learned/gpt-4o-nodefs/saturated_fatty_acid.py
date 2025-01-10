"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a long aliphatic chain with a carboxylic acid group at one end
    and no unsaturation or rings.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group at terminal position
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylic acid groups; requires exactly 1"
    
    # Check terminal carbon directly adjacent to carboxylate
    carboxylate_match = carboxylate_matches[0]
    carboxylic_carbon_idx = carboxylate_match[0]
    neighbors = mol.GetAtomWithIdx(carboxylic_carbon_idx).GetNeighbors()
    terminal = any(neighbor.GetDegree() == 1 for neighbor in neighbors)
    if not terminal:
        return False, "Carboxylic acid group is not at a terminal position"
    
    # Check for presence of any double bonds (unsaturation)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        return False, "Contains unsaturation (double bonds found)"
    
    # Check for rings
    if mol.HasSubstructMatch(Chem.MolFromSmarts("R")):
        return False, "Rings detected, should be aliphatic"
    
    # Verify chain length
    carbon_chain_pattern = Chem.MolFromSmarts("C" * 8)  # Minimum length of 8 for saturated fatty acids
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Chain length too short for fatty acid"

    # Ensure absence of other functional groups or heteroatoms besides oxygen
    non_carbon_heteroatoms = Chem.MolFromSmarts("[!C!O]")
    if mol.HasSubstructMatch(non_carbon_heteroatoms):
        return False, "Unexpected heteroatoms found"

    return True, "Molecule is a saturated fatty acid"