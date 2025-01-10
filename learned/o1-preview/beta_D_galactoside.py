"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is any D-galactoside having beta-configuration at its anomeric centre.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for beta-D-galactoside
    # This pattern represents the beta-D-galactopyranoside moiety with correct stereochemistry
    # The anomeric carbon (C1) is connected via an oxygen to any group ([*])
    beta_D_galactoside_smarts = "[C@H]1([O][*])O[C@H]([C@@H](O)[C@@H](O)[C@H]1O)CO"

    beta_D_galactoside_mol = Chem.MolFromSmarts(beta_D_galactoside_smarts)
    if beta_D_galactoside_mol is None:
        return False, "Failed to create beta-D-galactoside pattern"

    # Search for beta-D-galactoside substructure with chirality
    matches = mol.GetSubstructMatches(beta_D_galactoside_mol, useChirality=True)
    if matches:
        return True, "Contains beta-D-galactoside substructure"
    else:
        return False, "No beta-D-galactoside substructure found"