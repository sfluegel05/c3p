"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    Specifically, it checks for the presence of a beta-D-glucuronic acid moiety linked covalently via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for beta-D-glucuronic acid moiety with potential stereochemistry variations
    glucuronic_acid_pattern = Chem.MolFromSmarts("O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")

    # SMARTS pattern for glycosidic bond linkage, looking for any oxygen linkage
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@@H]1O[C@H](O)[C@@H](O)[C@H]1O[C@H]")

    # Check for the combined glucuronic acid structure and its glycosidic linkage correctly
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No compatible beta-D-glucuronic acid moiety found"

    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No proper glycosidic linkage found to beta-D-glucuronic acid"

    return True, "Contains beta-D-glucuronic acid moiety with correct glycosidic linkage"