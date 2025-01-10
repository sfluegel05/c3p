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
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (False, "Invalid SMILES string")

        # Define SMARTS pattern for beta-D-galactoside with beta configuration at anomeric carbon
        beta_d_galactoside_pattern = Chem.MolFromSmarts(
            "[C@H]1(O[C@H]([C@H]([C@H]([C@@H]([C@H]1O)O)O)O)CO)[O]*"
        )

        # Check if the molecule matches the beta-D-galactoside pattern
        if mol.HasSubstructMatch(beta_d_galactoside_pattern):
            return (True, "Contains beta-D-galactoside structural pattern with beta configuration")
        else:
            return (False, "No beta-D-galactoside pattern with correct beta-configuration found")
    
    except Exception as e:
        return (None, f"Error in processing: {str(e)}")

# Example usage
smiles_example = "OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)oc(-c2ccc(O)cc2)c(O)c3=O)[C@H](O)[C@@H](O)[C@H]1O"
print(is_beta_D_galactoside(smiles_example))