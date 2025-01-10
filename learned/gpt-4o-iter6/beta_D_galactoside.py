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

        # Define SMARTS pattern for beta-D-galactose ring with beta-linkage at anomeric center
        beta_galactose_pattern = Chem.MolFromSmarts('C1(O)C(O)C(O)C[C@@H](O)O1')
        
        # Check if the molecule contains the beta-D-galactose pattern
        if not mol.HasSubstructMatch(beta_galactose_pattern):
            return (False, "No beta-D-galactose structural pattern found")

        # Define SMARTS pattern to check for typical beta glycosidic linkage: attached -O-R
        glycosidic_pattern = Chem.MolFromSmarts('[C@H]1(O[*])')  # Looking for O-substituent at anomeric center
        if not mol.GetSubstructMatches(glycosidic_pattern):
            return (False, "Anomeric center lacks characteristic glycosidic bond")

        return (True, "Contains beta-D-galactoside structural pattern with beta configuration")

    except Exception as e:
        return (None, f"Error in processing: {str(e)}")

# Example usage
smiles_example = "OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)oc(-c2ccc(O)cc2)c(O)c3=O)[C@H](O)[C@@H](O)[C@H]1O"
print(is_beta_D_galactoside(smiles_example))