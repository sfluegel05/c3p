"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is defined as any D-galactoside having beta-configuration at its anomeric centre.

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

    # Identify the beta-D-galactoside core structure
    # The core can be approximated by the galactopyranose ring pattern
    beta_d_galactoside_pattern = Chem.MolFromSmarts(
        "[C@H]1([C@@H]([C@H](O)[C@@H]([C@@H]([C@H]1O)O)O)O)C"  # Represents beta-D-galactopyranoside
    )

    # Check if the molecule matches the beta-D-galactoside pattern
    if not mol.HasSubstructMatch(beta_d_galactoside_pattern):
        return False, "No beta-D-galactoside pattern found"
    
    return True, "Contains beta-D-galactoside pattern"

# Example usage
smiles_example = "OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)oc(-c2ccc(O)cc2)c(O)c3=O)[C@H](O)[C@@H](O)[C@H]1O"
print(is_beta_D_galactoside(smiles_example))