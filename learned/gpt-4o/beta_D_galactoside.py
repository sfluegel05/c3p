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
    
    # Identify the galactose moiety
    galactopyranoside_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@@H](O)[C@@H](O[C@H](CO)O1)O')
    
    if not mol.HasSubstructMatch(galactopyranoside_pattern):
        return False, "Galactose moiety not found"

    # Identify the beta anomeric configuration
    # The beta-configuration is given by [C@@H] next to the glycosidic oxygen and the anomeric carbon.
    beta_glycosidic_pattern = Chem.MolFromSmarts('[C@H]1(O[C@H](OC)CO)[C@@H](O)[C@@H](O)[C@@H]1O')
    
    if not mol.HasSubstructMatch(beta_glycosidic_pattern):
        return False, "Beta configuration at anomeric center not found"

    return True, "Beta-D-galactoside moiety found with proper beta configuration"