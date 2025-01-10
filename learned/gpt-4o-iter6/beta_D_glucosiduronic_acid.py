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

    # Broaden the search pattern for beta-D-glucuronic acid moiety (considering flexibility in stereochemistry and attaching points)
    # Central pattern in D-glucuronic acid with possible glycosidic and carboxy terminal
    glucuronic_acid_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H]([C@H](O)[C@@H](O)CO1)C(=O)O")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Broaden glycosidic linkage pattern to account for diversity in bond attachment
    # Checks if there is a heteroatom linkage to the sugar
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@H]1O[C@H](O)[C@@H](O)[C@H]1O[!H]")  # Make sure to capture ring-opening connections too
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No appropriate glycosidic linkage found"

    return True, "Contains the required beta-D-glucuronic acid moiety with correct linkage form"