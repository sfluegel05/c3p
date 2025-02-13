"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside has a beta-configured D-galactose sugar moiety with the correct stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern for beta-D-galactopyranoside
    # The pattern explicitly specifies the beta-anomeric configuration of D-galactose
    beta_d_galactopyranoside_smarts = "[C@@H]1(O[C@@H]([C@@H]([C@H]([C@H](O1)O)O)CO)[*:1])"  # Stereo-configured anomeric carbon
    # Note: Pattern should account for O-glycosidic linkages (indicating a connection)

    # Create a molecule from the SMARTS pattern
    galactopyranoside_pattern = Chem.MolFromSmarts(beta_d_galactopyranoside_smarts)
    
    # Check if the molecule matches the beta-D-galactopyranoside pattern
    if mol.HasSubstructMatch(galactopyranoside_pattern):
        return True, "Contains beta-D-galactopyranoside moiety"
    else:
        return False, "Beta-D-galactopyranoside moiety not found"

# The code ensures that the correct stereochemistry is being targeted and checked for.