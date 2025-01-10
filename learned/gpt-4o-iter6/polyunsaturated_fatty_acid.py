"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is identified as having a carboxyl group and more than one C=C bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for double bonds (C=C)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) <= 1:
        return False, f"Found {len(double_bond_matches)} double bonds, need more than 1"
    
    # Count carbon atoms to ensure it's a fatty acid (long chain of carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, f"Too few carbons ({c_count}) for a long chain fatty acid"
    
    return True, "Contains carboxylic acid group and more than one double bond, fitting the polyunsaturated fatty acid definition"