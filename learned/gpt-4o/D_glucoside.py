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

    # Define SMARTS patterns to recognize D-glucose structure in different stereochemical contexts
    # Including both alpha and beta anomers
    alpha_glucose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)[C@H]1O")
    beta_glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)[C@H]1O")

    # Check for glucose fingerprint in either anomeric form
    if not (mol.HasSubstructMatch(alpha_glucose_pattern) or mol.HasSubstructMatch(beta_glucose_pattern)):
        return False, "No D-glucose unit found"

    # Define a generic glycosidic linkage pattern; look for oxygen linkages involving sugar
    glycosidic_linkage_pattern = Chem.MolFromSmarts("C-O-[C@H]1O[C@@H](CO)[C@@H](O)[C@@H](O)[C@H]1O") 

    # Check for glycosidic linkage presence
    # Note: We look for broader connectivity that encompasses known D-glucose connection points
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "D-glucose unit found but no glycosidic linkage identified"

    return True, "Contains D-glucose unit with appropriate glycosidic linkage"