"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: CHEBI:35681 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid has two carboxylic acid groups and a hydroxy group
    on the alpha carbon of one of the carboxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if len(carboxyl_matches) != 2:
        return False, f"Found {len(carboxyl_matches)} carboxylic acid groups, need exactly 2"

    # Pattern for alpha-hydroxy carboxylic acid structure
    # Matches a carbon that has both an OH and a COOH group attached
    # Two possible arrangements are checked
    pattern1 = Chem.MolFromSmarts('[OH1]-[C](-[*,#1;!$(C(=O)[OH1])])-C(=O)[OH1]')  # OH-C-COOH
    pattern2 = Chem.MolFromSmarts('[OH1]-[C](-C(=O)[OH1])-[*,#1;!$(C(=O)[OH1])]')  # HOOC-C-OH
    
    # Check for either pattern
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    
    if not (matches1 or matches2):
        return False, "No alpha-hydroxy carboxylic acid group found"
    
    # Verify that the second carboxylic acid group is not directly attached to the alpha carbon
    alpha_carbons = set()
    for match in matches1:
        alpha_carbons.add(match[1])  # Index 1 is the alpha carbon
    for match in matches2:
        alpha_carbons.add(match[1])
        
    # For each potential alpha carbon, verify it has exactly one COOH group
    for alpha_carbon_idx in alpha_carbons:
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        cooh_count = 0
        
        # Count COOH groups attached to this carbon
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                # Check if this carbon is part of a COOH group
                if any(alpha_carbon_idx == match[1] and neighbor.GetIdx() == match[2] 
                      for match in carboxyl_matches):
                    cooh_count += 1
        
        if cooh_count == 1:
            return True, "Contains two carboxylic acid groups with a hydroxy group on the alpha carbon"
            
    return False, "No valid alpha-hydroxy carboxylic acid structure found"