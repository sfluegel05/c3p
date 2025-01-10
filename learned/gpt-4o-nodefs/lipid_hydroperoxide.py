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
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a hydroperoxide group (-O-O-H)
    hydroperoxide_pattern = Chem.MolFromSmarts("O[OH]")
    if not mol.HasSubstructMatch(hydroperoxide_pattern):
        return False, "No hydroperoxide group found"
    
    # Check for the presence of the carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbon atoms to assess if the molecule is of a lipid-like size (≥14 carbons for revisiting)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:  # Adjusted to accommodate slightly shorter lipids
        return False, "Carbon chain is too short for a lipid-like structure"
    
    # Check for multiple double bonds to ensure polyunsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bonds < 2:  # More specific criterion for polyunsaturation
        return False, f"Insufficient double bonds for polyunsaturation, found {double_bonds}"
        
    return True, "Molecule matches structure of a lipid hydroperoxide"