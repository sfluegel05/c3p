"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside contains a D-glucose unit which can be identified with generalized SMARTS patterns.

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

    # Define SMARTS patterns to recognize D-glucose or its derivatives
    glucose_patterns = [
        # D-glucose or variations in open/closed form contexts
        Chem.MolFromSmarts("C1([C@@H]([C@H]([C@H]([C@@H](C1O)O)O)O)O)"),
        Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H](C(O)O)[C@H](O1)CO)]"),
        # Including terms for derivatives (e.g., glucopyranoside, glucosidic forms)
        Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)[C@H]1O"),
        Chem.MolFromSmarts("C[C@@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)[C@H]1O"),
        # Glucuronic acid version: carboxyl group included
        Chem.MolFromSmarts("OC(=O)[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"),
    ]

    # Check for the presence of any defined glucose patterns
    glucose_found = False
    for pattern in glucose_patterns:
        if mol.HasSubstructMatch(pattern):
            glucose_found = True
            break

    if not glucose_found:
        return False, "No D-glucose unit or derivative found"

    # Define more general glycosidic linkage pattern
    general_glycosidic_linkage = Chem.MolFromSmarts("C-O[C@H]")

    # Check for any acceptable glycosidic linkage involving a glucose derivative
    if not mol.HasSubstructMatch(general_glycosidic_linkage):
        return False, "D-glucose derivative found but no glycosidic linkage identified"

    return True, "Contains D-glucose unit or derivative with appropriate glycosidic linkage"