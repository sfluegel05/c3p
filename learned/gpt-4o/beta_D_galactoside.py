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
    
    # More flexible galactose moiety pattern to detect the presence of galactose
    galactose_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O')
    
    # Identify the galactose moiety
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "Galactose moiety not found"

    # Identify the beta configuration at the anomeric center
    anomeric_c = Chem.MolFromSmarts('[C@H]1(O[C@H](CO)*)O[C@H](O)[C@@H](O)[C@H]1O')
    beta_match = mol.HasSubstructMatch(anomeric_c)

    if not beta_match:
        return False, "Beta configuration at anomeric center not found"

    return True, "Beta-D-galactoside moiety found with proper beta configuration"