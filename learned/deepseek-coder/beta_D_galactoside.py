"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: CHEBI:18133 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is a galactose ring in the beta-configuration, where the hydroxyl group at the anomeric carbon (C1) is in the axial position.

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

    # Define the beta-D-galactoside pattern
    # The pattern matches a galactose ring with the beta-configuration at the anomeric center
    beta_D_galactoside_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    
    # Check if the molecule contains the beta-D-galactoside pattern
    if mol.HasSubstructMatch(beta_D_galactoside_pattern):
        return True, "Contains a galactose ring in the beta-configuration"
    else:
        return False, "No galactose ring in the beta-configuration found"

# Example usage:
# smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"  # Methyl beta-D-galactoside
# print(is_beta_D_galactoside(smiles))