"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-glucoside
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is any glucoside in which the glycoside group is derived from D-glucose.

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

    # Define SMARTS pattern for D-glucose moiety with glycosidic linkage
    # at the anomeric carbon (C1). The '*' allows for any attachment.
    d_glucoside_smarts = "[C@@H]1(O[*])[C@H](O)[C@@H](O)[C@H](O)[C@H]1CO"

    d_glucoside_mol = Chem.MolFromSmarts(d_glucoside_smarts)
    if d_glucoside_mol is None:
        return False, "Failed to generate D-glucoside SMARTS pattern"

    # Use substructure search to find D-glucose moiety with glycosidic linkage
    matches = mol.GetSubstructMatches(d_glucoside_mol, useChirality=True)
    if not matches:
        return False, "No D-glucoside moiety found"

    return True, "Contains D-glucose moiety connected via glycosidic bond"