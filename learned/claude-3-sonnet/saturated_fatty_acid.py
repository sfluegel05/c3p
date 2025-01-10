"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: saturated fatty acids
Definition: Any fatty acid containing no carbon to carbon multiple bonds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should only have one
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Check for absence of carbon-carbon multiple bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    
    if mol.HasSubstructMatch(double_bond_pattern):
        return False, "Contains carbon-carbon double bonds"
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon-carbon triple bonds"
    
    # Count carbons and check chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 2:
        return False, "Carbon chain too short for fatty acid"
        
    # Check atom types (allowing for deuterium)
    allowed_atoms = {1, 6, 8}  # H, C, O
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains atoms other than C, H, O"
        
    # Check that carboxylic acid is terminal
    # First get the carbon of the carboxyl group
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_matches[0][0])
    
    # Count non-oxygen connections to carboxyl carbon
    non_oxygen_connections = sum(1 for neighbor in carboxyl_carbon.GetNeighbors() 
                               if neighbor.GetAtomicNum() != 8)
    
    if non_oxygen_connections > 1:
        return False, "Carboxylic acid group is not terminal"
        
    # Count oxygens - should have exactly 2 for the carboxyl group
    # (plus any additional OH groups in hydroxy fatty acids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for carboxylic acid"
        
    return True, "Saturated fatty acid with terminal carboxyl group"