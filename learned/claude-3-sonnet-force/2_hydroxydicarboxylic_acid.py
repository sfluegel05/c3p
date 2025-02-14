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
    
    # Check for two carboxylic acid groups (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[$(C(=O)O)]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 2:
        return False, "Not a dicarboxylic acid"
    
    # Check for hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"
    
    # Find alpha carbon (connected to hydroxyl group)
    alpha_carbon_idx = None
    for match in hydroxyl_matches:
        hydroxyl_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in hydroxyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbon_idx = neighbor.GetIdx()
                break
        if alpha_carbon_idx is not None:
            break
    
    if alpha_carbon_idx is None:
        return False, "No alpha carbon found connected to hydroxyl group"
    
    # Check if alpha carbon is connected to one of the carboxylic acid groups
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    for neighbor in alpha_carbon.GetNeighbors():
        for bond in neighbor.GetBonds():
            if bond.GetBeginAtomIdx() == neighbor.GetIdx():
                continue
            if mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomicNum() == 6:
                if mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomicNum() == 8:
                    if mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetTotalNumHs() == 0:
                        return True, "2-hydroxydicarboxylic acid"
    
    return False, "Alpha carbon not connected to carboxylic acid group"