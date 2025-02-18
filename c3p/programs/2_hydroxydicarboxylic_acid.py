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
    carboxyl_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if len(carboxyl_matches) != 2:
        return False, f"Found {len(carboxyl_matches)} carboxylic acid groups, need exactly 2"

    # Find alpha-hydroxy carboxylic acid pattern
    # [OH]-[C](-[*])(-[*])-C(=O)-[OH]
    alpha_hydroxy_acid_pattern = Chem.MolFromSmarts('[OX2H1][CX4]([*])([*])C(=[OX1])[OX2H1]')
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_acid_pattern)
    
    if not alpha_hydroxy_matches:
        return False, "No alpha-hydroxy carboxylic acid group found"

    # Additional check to ensure the hydroxy group is actually on an alpha carbon
    # of one of the carboxyl groups (not just anywhere in the molecule)
    alpha_carbon_atoms = set()
    for match in carboxyl_matches:
        carboxyl_carbon = match[0]  # First atom in carboxyl match is the carbon
        for neighbor in mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                alpha_carbon_atoms.add(neighbor.GetIdx())

    hydroxy_carbons = set()
    for match in alpha_hydroxy_matches:
        hydroxy_carbons.add(match[1])  # Second atom in match is the carbon with OH

    if not (alpha_carbon_atoms & hydroxy_carbons):
        return False, "Hydroxy group not on alpha carbon of carboxyl group"

    return True, "Contains two carboxylic acid groups with a hydroxy group on the alpha carbon"