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

    # Look for beta-D-glucuronic acid moiety pattern
    glucuronic_pattern = Chem.MolFromSmarts("[C@H]1OC[C@@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glucuronic_pattern):
        return False, "No beta-D-glucuronic acid moiety found"
    
    # Look for carboxylic acid group [C(=O)O]
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group attached to the glucuronic acid moiety"
 
    # Check for glycosidic linkage -O- connectivity from glucuronic acid
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H]1OC[C@@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Missing glycosidic bond attached to glucuronic acid"

    # Look for anomeric carbon linkage via glycosidic bond
    anomeric_link_pattern = Chem.MolFromSmarts("[C@H]1OC[C@@H](O)[C@@H](O)[C@H](O1)-O")
    if not mol.HasSubstructMatch(anomeric_link_pattern):
        return False, "Anomeric carbon linkage via a glycosidic bond not found"

    return True, "Contains beta-D-glucuronic acid moiety attached via glycosidic bond"