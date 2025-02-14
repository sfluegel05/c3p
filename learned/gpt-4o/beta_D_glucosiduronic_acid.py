"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid involves a beta-D-glucuronic acid connected glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a molecule is beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Define a more general SMARTS for glucuronic acid backbone
    glucuronic_acid_smarts = "OC[C@@H]1OC(O)[C@H](O[C@H]1O)C(=O)O"
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)

    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid substructure found"

    # Flexible pattern for glycosidic linkage
    glycosidic_linkage_smarts = "O[C@H]"
    glycosidic_pattern = Chem.MolFromSmarts(glycosidic_linkage_smarts)

    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"
    
    return True, "Detected beta-D-glucuronic acid with a glycosidic linkage"

# Example usage:
# is_beta_D_glucosiduronic_acid("...")