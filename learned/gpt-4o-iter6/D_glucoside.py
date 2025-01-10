"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside has a D-glucose moiety with a glycosidic linkage to another molecule.
    
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

    # Define SMARTS pattern for the beta-D-glucose moiety
    d_glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    # Look for glycosidic linkage, typically an ether-like linkage from the anomeric carbon
    glycosidic_oxygen_pattern = Chem.MolFromSmarts("CO[C@@H]1O[C@H]")  # Adjust pattern for better coverage

    # Check for D-glucose moiety with potential glycosidic linkage
    if mol.HasSubstructMatch(d_glucose_pattern) and mol.HasSubstructMatch(glycosidic_oxygen_pattern):
        return True, "Contains D-glucopyranosyl moiety with potential glycosidic linkage"

    return False, "Does not contain the D-glucopyranosyl moiety with required stereochemistry or linkage"