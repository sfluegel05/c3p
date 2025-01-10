"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a precise hydroperoxide group pattern
    hydroperoxide_pattern = Chem.MolFromSmarts("O[OH]")
    if not mol.HasSubstructMatch(hydroperoxide_pattern):
        return False, "No hydroperoxide group found"
    
    # Define a carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Ensure long carbon chain presence (minimum 12 carbons for lipid-like)
    carbon_chain = Chem.MolFromSmarts("[C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Carbon chain is too short for a lipid"
    
    # Check for sufficient polyunsaturation (at least 2 double bonds)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bonds) < 2:
        return False, "Insufficient double bonds"
        
    return True, "Molecule matches structure of a lipid hydroperoxide"