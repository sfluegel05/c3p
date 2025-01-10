"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a glucoside in which the glycoside group is derived from D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # D-glucose: a pyranose ring with specific hydroxyl stereochemistry
    d_glucose_pattern1 = Chem.MolFromSmarts('OC[C@H]1O[C@@H]([C@H]([C@@H](O)[C@@H]1O)O)CO')  # Beta-D-glucose
    d_glucose_pattern2 = Chem.MolFromSmarts('OC[C@@H]1O[C@H]([C@@H]([C@H](O)[C@H]1O)O)CO')  # Alpha-D-glucose

    # Check for D-glucose moiety
    if not (mol.HasSubstructMatch(d_glucose_pattern1) or mol.HasSubstructMatch(d_glucose_pattern2)):
        return False, "D-glucose moiety not found"

    # Check for glycosidic bonds
    # Looking for any glycosidic linkage through the anomeric carbon (O-C1-O)
    glycosidic_bond_pattern1 = Chem.MolFromSmarts('[C@H]1(O)[C@@H]([C@H](O)O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)OC(O)[C@H]1O')
    glycosidic_bond_pattern2 = Chem.MolFromSmarts('[C@@H]1(O)[C@H]([C@@H](O)O[C@H]2O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]2O)OC(O)[C@@H]1O')

    if mol.HasSubstructMatch(glycosidic_bond_pattern1) or mol.HasSubstructMatch(glycosidic_bond_pattern2):
        return True, "Contains D-glucose moiety with glycosidic linkage"

    return False, "Glycosidic linkage not found to D-glucose moiety"