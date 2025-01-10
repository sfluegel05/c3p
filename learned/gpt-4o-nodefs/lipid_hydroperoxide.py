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
    
    # Look for the hydroperoxide group (-O-O-H)
    hydroperoxide_pattern = Chem.MolFromSmarts("O[OH]")
    if not mol.HasSubstructMatch(hydroperoxide_pattern):
        return False, "No hydroperoxide group found"
    
    # Look for the carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count the total number of carbon atoms to check for lipid-like size (≥16 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:
        return False, "Carbon chain is too short for a lipid-like structure"
    
    # Check for polyunsaturation (presence of any double bond)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bonds found for polyunsaturation"
        
    return True, "Molecule matches structure of a lipid hydroperoxide"