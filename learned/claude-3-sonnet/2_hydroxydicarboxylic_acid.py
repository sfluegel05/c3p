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

    # Pattern for alpha-hydroxy carboxylic acid:
    # [OH]-[C]-C(=O)[OH] where the central carbon is connected to exactly one OH
    # and one COOH, ensuring it's truly the alpha position
    alpha_hydroxy_pattern = Chem.MolFromSmarts('[OH1]-[C;!$(C-[OH1])](-[*,#1])-C(=O)[OH1]')
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    
    if not alpha_hydroxy_matches:
        return False, "No alpha-hydroxy carboxylic acid group found"

    # For each alpha-hydroxy match, verify it's a true alpha position
    valid_alpha = False
    for match in alpha_hydroxy_matches:
        alpha_carbon = mol.GetAtomWithIdx(match[1])  # The central carbon atom
        
        # Count number of COOH groups attached to this carbon
        cooh_count = 0
        oh_count = 0
        for neighbor in alpha_carbon.GetNeighbors():
            # Check if neighbor is part of COOH
            if neighbor.GetSymbol() == 'C':
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'O' and nn.GetTotalNumHs() == 1:
                        cooh_count += 1
            # Check if neighbor is OH
            elif neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                oh_count += 1
                
        if cooh_count == 1 and oh_count == 1:
            valid_alpha = True
            break

    if not valid_alpha:
        return False, "No valid alpha-hydroxy carboxylic acid group found"

    return True, "Contains two carboxylic acid groups with a hydroxy group on the alpha carbon"