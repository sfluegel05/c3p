"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-D-galactoside unit
    beta_D_galactoside_smarts = "[C@H]1([O][*])[O][C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO"
    beta_D_galactoside_pattern = Chem.MolFromSmarts(beta_D_galactoside_smarts)
    if beta_D_galactoside_pattern is None:
        return False, "Failed to create beta-D-galactoside pattern"

    # Search for beta-D-galactoside substructure
    if mol.HasSubstructMatch(beta_D_galactoside_pattern):
        return True, "Contains beta-D-galactoside substructure"
    else:
        return False, "No beta-D-galactoside substructure found"