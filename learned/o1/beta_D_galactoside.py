"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is a D-galactose residue with beta-configuration at its anomeric center
    connected via a glycosidic bond to any aglycone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for beta-D-galactose with a glycosidic linkage
    # Stereochemistry for D-galactose:
    # C1 (anomeric carbon): beta (up)
    # C2: R configuration
    # C3: S configuration
    # C4: R configuration
    # C5: R configuration
    # The anomeric oxygen is connected to an aglycone (ANY group except hydrogen)
    beta_D_galactoside_smarts = """
    [C@H]1([O;X2;$([OH]),$([O][#6]),$([O][#7]),$([O][#8]),$([O][#15])])  # C1 with beta configuration and glycosidic oxygen
    [C@@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O                           # Ring atoms with specific stereochemistry
    """

    # Remove whitespace and newlines
    beta_D_galactoside_smarts = ''.join(beta_D_galactoside_smarts.split())

    beta_D_galactoside_mol = Chem.MolFromSmarts(beta_D_galactoside_smarts)
    if beta_D_galactoside_mol is None:
        return False, "Error in SMARTS pattern for beta-D-galactoside"

    # Check if the molecule contains the beta-D-galactoside substructure
    matches = mol.GetSubstructMatches(beta_D_galactoside_mol, useChirality=True)
    if matches:
        # Further verify that the anomeric oxygen is part of a glycosidic bond
        for match in matches:
            anomeric_carbon_idx = match[0]
            anomeric_oxygen_idx = mol.GetAtomWithIdx(anomeric_carbon_idx).GetNeighbors()[0].GetIdx()
            # Ensure anomeric oxygen is not bound to hydrogen (not an OH group)
            anomeric_oxygen = mol.GetAtomWithIdx(anomeric_oxygen_idx)
            if anomeric_oxygen.GetAtomicNum() == 8 and anomeric_oxygen.GetTotalNumHs() == 0:
                return True, "Contains beta-D-galactoside substructure"
        return False, "Anomeric oxygen is not part of a glycosidic bond"
    else:
        return False, "Does not contain beta-D-galactoside substructure"