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

    # Define SMARTS pattern for the D-glucose moiety with both alpha and beta forms
    d_glucose_pattern_beta = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1CO")
    d_glucose_pattern_alpha = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1CO")
    glycosidic_oxygen = Chem.MolFromSmarts("CO[C@@H]1O")

    # Check for D-glucose moiety with possible glycosidic linkage
    if (mol.HasSubstructMatch(d_glucose_pattern_beta) or mol.HasSubstructMatch(d_glucose_pattern_alpha)) and mol.HasSubstructMatch(glycosidic_oxygen):
        return True, "Contains D-glucopyranosyl moiety with potential glycosidic linkage"

    return False, "Does not contain the D-glucopyranosyl moiety with required stereochemistry or linkage"