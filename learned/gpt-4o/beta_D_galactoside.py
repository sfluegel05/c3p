"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A D-galactoside has the galactose moiety with a beta-configuration at its anomeric center.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated and more flexible pattern to detect a beta-D-galactose moiety
    galactose_pattern = Chem.MolFromSmarts('[C@H]1([O][C@@H]([C@H](O)[C@H]([C@H]1O)O)O)CO')
    
    # Match the galactose moiety with beta configuration
    matches = mol.GetSubstructMatches(galactose_pattern)

    if not matches:
        return False, "Beta configuration at anomeric center not found or galactose moiety not identified"
    
    return True, "Beta-D-galactoside moiety found with proper beta configuration"