"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for beta-hydroxy pattern (OH on carbon 3 positions away from COOH)
    # Pattern: HO-C-C-C(=O)OH
    beta_hydroxy_pattern = Chem.MolFromSmarts("[OX2H1]-[CX4]-[CX4]-[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No hydroxy group at beta position (3-position)"

    # Count carbons to ensure it's a fatty acid (minimum 5 carbons for fatty acids)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, "Carbon chain too short to be a fatty acid"
    
    # Count rotatable bonds to verify chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too rigid to be a fatty acid"
    
    # Additional checks to exclude non-fatty acid structures
    # Check if molecule has reasonable number of oxygens (at least 3: COOH + OH)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms for 3-hydroxy fatty acid"
    
    # Verify the carbon chain is primarily linear (not too branched)
    # Count carbons with more than 2 carbon neighbors
    highly_branched_carbons = sum(1 for atom in mol.GetAtoms() 
                                if atom.GetAtomicNum() == 6 and 
                                len([n for n in atom.GetNeighbors() 
                                    if n.GetAtomicNum() == 6]) > 2)
    if highly_branched_carbons > carbon_count / 4:  # Allow some branching but not too much
        return False, "Carbon chain too branched for typical fatty acid"

    return True, "Contains carboxylic acid group and hydroxy group at 3-position with appropriate fatty acid chain"