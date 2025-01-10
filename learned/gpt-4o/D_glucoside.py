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
    
    # D-glucose substructure (common cyclic form with oxygens and hydroxyls)
    # D-glucose has a six-membered ring and hydroxyl group pattern
    d_glucose_pattern = Chem.MolFromSmarts('C1([C@@H]2[C@@H]([C@H](C([C@H](O1)O)O)O)O2)O') # Simplified representation

    # Check for D-glucose substructure
    if not mol.HasSubstructMatch(d_glucose_pattern):
        return False, "D-glucose moiety not found"

    # Glycosidic linkage pattern (oxygen connected to glucose ring and another moiety)
    glycosidic_linkage_pattern = Chem.MolFromSmarts('[O][C@H]1[C@@H]([C@H]([C@@H](O[C@@H]2O[C@H](O)[C@H](O)[C@@H](O)C2)O)O)O1')

    if mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return True, "Contains D-glucose moiety with glycosidic linkage"

    return False, "Glycosidic linkage not found to D-glucose moiety"