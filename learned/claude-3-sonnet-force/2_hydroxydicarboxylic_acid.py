"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: CHEBI:51139 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is any dicarboxylic acid carrying a hydroxy group on the
    carbon atom at position alpha to the carboxy group.

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
    carboxylic_acid_pattern = Chem.MolFromSmarts("[$(C(=O)O)]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 2:
        return False, "Not a dicarboxylic acid"
    
    # Find hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"
    
    # Check if any hydroxyl group is alpha to a carboxylic acid group
    for hydroxyl_idx in hydroxyl_matches:
        for acid_idx in carboxylic_acid_matches:
            path = Chem.FindShortestPath(mol, hydroxyl_idx[0], acid_idx[0])
            if len(path) == 3:  # Alpha carbon is in the middle
                alpha_carbon_idx = path[1]
                return True, "2-hydroxydicarboxylic acid"
    
    return False, "No hydroxyl group found alpha to a carboxylic acid group"