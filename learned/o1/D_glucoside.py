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

    # Define a SMARTS pattern for D-glucose moiety in pyranose form
    # This pattern considers both alpha and beta anomers and correct stereochemistry
    d_glucose_smarts = """
    [C@H]1([O])[C@@H]([C@H](O)[C@@H](O)[C@H](O)[C@H]1O)CO
    """

    d_glucose_mol = Chem.MolFromSmarts(d_glucose_smarts)
    if d_glucose_mol is None:
        return False, "Failed to generate D-glucose SMARTS pattern"

    # Use substructure search to find D-glucose moiety
    matches = mol.GetSubstructMatches(d_glucose_mol, useChirality=True)
    if not matches:
        return False, "No D-glucose moiety found"

    # For each match, check for glycosidic linkage at the anomeric carbon (C1)
    for match in matches:
        glucose_atoms = set(match)
        anomeric_carbon_idx = match[0]  # C1 in our SMARTS pattern
        anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)

        # Check bonds from the anomeric carbon
        glycosidic_bond_found = False
        for bond in anomeric_carbon.GetBonds():
            neighbor = bond.GetOtherAtom(anomeric_carbon)
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in glucose_atoms:
                # Atom outside the glucose moiety
                if neighbor.GetAtomicNum() in [8, 7]:  # Oxygen or Nitrogen
                    glycosidic_bond_found = True
                    break  # Glycosidic bond found
        if glycosidic_bond_found:
            return True, "Contains D-glucose moiety connected via glycosidic bond"

    return False, "No glycosidic linkage found at anomeric carbon"