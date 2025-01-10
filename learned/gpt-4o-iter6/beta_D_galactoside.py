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

        # Define SMARTS pattern for beta-D-galactose ring
        # Pattern includes: six-membered pyranose ring, hydroxyl groups in correct orientation
        beta_glucose_ring = Chem.MolFromSmarts('OC[C@@H]1O[C@H]([C@H](O)[C@@H](O)[C@H]1O)')

        # Check if the molecule contains the beta-D-galactose ring pattern
        if not mol.HasSubstructMatch(beta_glucose_ring):
            return (False, "No beta-D-galactose ring pattern found")

        # Extract potential glycosidic bonds (e.g., -O-R) that might be attached to the anomeric center
        anomeric_oxygen_pattern = Chem.MolFromSmarts("OC")

        if not mol.HasSubstructMatch(anomeric_oxygen_pattern):
            return (False, "No potential glycosidic bond at anomeric center")

        return (True, "Contains beta-D-galactoside structural pattern with beta configuration")
    
    except Exception as e:
        return (None, f"Error in processing: {str(e)}")

# Example usage
smiles_example = "OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)oc(-c2ccc(O)cc2)c(O)c3=O)[C@H](O)[C@@H](O)[C@H]1O"
print(is_beta_D_galactoside(smiles_example))