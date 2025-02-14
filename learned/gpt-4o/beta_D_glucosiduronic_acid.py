"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid involves a beta-D-glucuronic acid glycosidically linked to another moiety.

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

    # Define SMARTS patterns
    # Pattern for glucuronic acid with appropriate beta linkage
    glucuronic_acid_smarts = "[C@H]1(O)[C@@H](O)[C@H](O)[C@H](O1)C(=O)O"  # Stereochemistry judged for beta-D
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid substructure found"

    # Check for glycosidic linkage specifically in the beta position
    # We need to ensure it's a bond involving the non-reducing end of glucuronic acid
    beta_linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)")
    if not mol.HasSubstructMatch(beta_linkage_pattern):
        return False, "No beta glycosidic linkage found"
    
    return True, "Contains beta-D-glucuronic acid with glycosidic linkage"

# Example usage:
# is_beta_D_glucosiduronic_acid("...")