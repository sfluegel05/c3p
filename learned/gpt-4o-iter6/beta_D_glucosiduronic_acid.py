"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid has a beta-D-glucuronic acid moiety bound via a glycosidic bond.

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

    # Pattern for a glucuronic acid moiety allowing for flexibility in stereochemistry
    glucuronic_acid_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)CO1")
    
    # Combined pattern to check for the core glucuronic acid ring with the glycosidic bond
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@@H]1O[C@H](O)[C@H](O)[C@H]1O[C@H]")

    # Ensure the glucuronic acid core and its glycosidic connections are present
    if not mol.HasSubstructMatch(glucuronic_acid_pattern) or not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No beta-D-glucuronic acid moiety found or improper linkage"

    return True, "Contains the required beta-D-glucuronic acid moiety with correct linkage form"