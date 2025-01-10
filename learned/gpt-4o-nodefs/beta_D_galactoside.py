"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside contains the beta-D-galactopyranosyl unit characterized
    by specific stereochemistry around the pyranose ring.

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

    # Define an improved beta-D-galactopyranosyl pattern
    # Updated the stereochemistry pattern and added flexibility
    # to capture variations in substitution.
    beta_d_galactoside_pattern = Chem.MolFromSmarts(
        "[C@H]1(O)[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O1"
    )
    
    if mol.HasSubstructMatch(beta_d_galactoside_pattern):
        return True, "Contains beta-D-galactopyranosyl unit"

    return False, "Does not contain beta-D-galactopyranosyl unit"