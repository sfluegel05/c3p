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
    A straight-chain saturated fatty acid is defined as any saturated fatty acid lacking a side-chain.
    
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

    # Check for rings
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains cyclic structures"

    # Check for hydroxyl groups (except the one in carboxylic acid)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1][CX4]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Contains hydroxyl group(s)"

    # Check for ketone groups
    ketone_pattern = Chem.MolFromSmarts("[CX4]-[CX3](=O)-[CX4]")
    if mol.HasSubstructMatch(ketone_pattern):
        return False, "Contains ketone group(s)"
        
    # Get the carboxylic acid carbon
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_matches[0][0])
    
    # Traverse the carbon chain from the carboxylic group
    visited = set()
    current = carboxyl_carbon
    chain_length = 0
    branching_detected = False
    
    while current is not None:
        visited.add(current.GetIdx())
        chain_length += 1
        
        # Find next carbon in chain
        next_carbon = None
        carbon_neighbors = 0
        
        for neighbor in current.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                carbon_neighbors += 1
                next_carbon = neighbor
                
        # Check for branching
        if carbon_neighbors > 1:
            branching_detected = True
            break
            
        current = next_carbon

    if branching_detected:
        return False, "Branched carbon chain detected (side chain present)"

    # Check minimum chain length (typically 4+ carbons for fatty acids)
    if chain_length < 4:
        return False, "Carbon chain too short for a fatty acid"

    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(double_bond_pattern):
        return False, "Contains double bonds (unsaturated)"
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains triple bonds (unsaturated)"

    # Check that all carbons in the main chain are saturated (except carboxylic carbon)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Skip the carboxylic carbon
            if atom.GetIdx() == carboxyl_matches[0][0]:
                continue
            # Check that all other carbons are sp3 hybridized
            if atom.GetHybridization() != Chem.HybridizationType.SP3:
                return False, "Contains unsaturated carbons"

    # All checks passed
    return True, "Straight-chain saturated fatty acid with single carboxylic acid group"