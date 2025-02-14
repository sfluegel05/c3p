"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid with a hydroxy group on the alpha carbon to one of the carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    hydroxy_pattern = Chem.MolFromSmarts("CO")  # Generic C-OH pattern
    carbon_pattern = Chem.MolFromSmarts("[CX4;H2,H1]")  # Sp3 carbon

    # Find all carboxyl groups
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)    
    if len(carboxylic_matches) < 2:
        return False, f"Found {len(carboxylic_matches)} carboxylic acid groups, need exactly 2"

    # Check for alpha-hydroxy relation
    for match in carboxylic_matches:
        carboxy_carbon_idx = match[0]

        # Check for an alpha-hydroxy group
        for neighbor in mol.GetAtomWithIdx(carboxy_carbon_idx).GetNeighbors():
            if neighbor.GetIdx() == match[1]:  # Skip the carbonyl oxygen attached to the carbon
                continue
            
            if mol.HasSubstructMatch(hydroxy_pattern, useChirality=True, atomMap={0: neighbor.GetIdx()}):
                return True, "Contains two carboxylic groups with a hydroxy group on alpha carbon"

    return False, "Does not have a hydroxy group on an alpha carbon to a carboxylic acid group"