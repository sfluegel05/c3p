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
    
    # Refined galactose moiety pattern
    galactopyranoside_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@@H](O)[C@@H](O[C@@H]([C@H]([C@@H]1O)O)CO)O')
    
    # Identify the galactose moiety
    if not mol.HasSubstructMatch(galactopyranoside_pattern):
        return False, "Galactose moiety not found"

    # Improved beta anomeric configuration pattern
    beta_glycosidic_pattern = Chem.MolFromSmarts('[C@@H]1(O[C@@H](OC)CO)[C@H](O)[C@H](O)[C@@H]1O') 
    
    # Identify the beta anomeric configuration
    if not mol.HasSubstructMatch(beta_glycosidic_pattern):
        return False, "Beta configuration at anomeric center not found"

    return True, "Beta-D-galactoside moiety found with proper beta configuration"