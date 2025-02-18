"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: straight-chain saturated fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should only have one
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
        
    # Check for any branching (carbon atoms with more than 2 carbon neighbors)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            carbon_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
            if carbon_neighbors > 2:
                return False, "Branched carbon chain detected"
    
    # Check for unsaturation (double or triple bonds)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    
    if mol.HasSubstructMatch(double_bond_pattern):
        return False, "Contains double bonds (unsaturated)"
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains triple bonds (unsaturated)"
    
    # Count carbons in the main chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 2:
        return False, "Carbon chain too short for a fatty acid"
        
    # Check for non C,H,O atoms (excluding deuterium)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:  # H, C, O
            if not (atomic_num == 1 and atom.GetIsotope() == 2):  # Allow deuterium
                return False, "Contains elements other than C, H, O (or deuterium)"
    
    # Check for hydroxyl groups on the chain (excluding the carboxylic acid)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Contains hydroxyl groups on the carbon chain"
    
    # Check for ketone groups
    ketone_pattern = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4]")
    if mol.HasSubstructMatch(ketone_pattern):
        return False, "Contains ketone groups"
    
    # All checks passed
    return True, "Straight-chain saturated fatty acid with single carboxylic acid group"