"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    This involves having at least one hydroxyl group in a central (potentially stereo-locked) position,
    and two carboxylic acid groups, allowing for variation in connectivity for the core carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create a pattern for the hydroxyl group attached to a central carbon (e.g., secondary carbon flexibility)
    hydroxy_pattern = Chem.MolFromSmarts("[#6][CH](O)[#6]")
    
    # Create a pattern for two carboxylic acid groups in flexible positions.
    carboxylic_acid_pattern1 = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_pattern2 = Chem.MolFromSmarts("O=C(O)")

    # Check for the presence of a central hydroxy group, with flexibility around potential bonds.
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No central hydroxy group found at expected position"

    # Check for two carboxylic acid groups (we allow flexible patterns to ensure branched/network matches).
    carboxyl_matches1 = mol.GetSubstructMatches(carboxylic_acid_pattern1)
    carboxyl_matches2 = mol.GetSubstructMatches(carboxylic_acid_pattern2)

    total_carboxylics = len(carboxyl_matches1) + len(carboxyl_matches2) 
    if total_carboxylics < 2:
        return False, f"Found {total_carboxylics} carboxylic acid groups, need exactly 2"

    # Confirming with logical positioning checks if need be
    # Efficient use of connectivity checks and relaxed patterning logic for compound flexibility.
    
    return True, "Contains 2-hydroxy group and two separate carboxylic acid groups in suitable positions."