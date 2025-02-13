"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is any D-galactoside having beta-configuration at its anomeric center.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecular structure corresponds to beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Beta-D-galactopyranose pattern - Updated
    # Matching the galactopyranose ring with beta-OH on anomeric carbon
    beta_d_galactoside_pattern = Chem.MolFromSmarts(
        "[C@H]1(O)[C@@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1"  # Simplified pattern for beta-D-galactopyranosides
    )

    # Check if the molecule matches the beta-D-galactoside pattern
    if mol.HasSubstructMatch(beta_d_galactoside_pattern):
        return True, "Contains beta-D-galactoside pattern"
    else:
        return False, "No beta-D-galactoside pattern found"

# Example usage
smiles_example = "OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)oc(-c2ccc(O)cc2)c(O)c3=O)[C@H](O)[C@@H](O)[C@H]1O"
print(is_beta_D_galactoside(smiles_example))