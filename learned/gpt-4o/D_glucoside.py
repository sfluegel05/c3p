"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside contains a D-glucose unit attached via a glycoside bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for D-glucose based unit, covering both alpha and beta anomers
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)[C@H]1O")
    
    # Check if the molecule contains at least one D-glucose unit
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No D-glucose unit found"
    
    # Define a more generalized pattern for glycosidic linkage
    # This pattern attempts to match an O-glycosidic linkage with at least one D-glucose unit
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[C,c,O]O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)[C@H]1O")
    
    # Check if the molecule contains a glycosidic linkage
    # Note: This is a generalized linkage pattern
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage found"

    return True, "Contains D-glucose unit with glycosidic linkage"