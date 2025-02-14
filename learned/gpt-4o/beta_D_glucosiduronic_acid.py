"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid involves a beta-D-glucuronic acid connected glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a molecule is beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general SMARTS pattern for glucuronic acid backbone
    # Pattern: C1(OC(C)[O])OC(O)[C@@H](O)C1C(=O)O
    glucuronic_acid_smarts = "OC1C(O)C(O[C@H]1O)COC(=O)[CH2]"
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)

    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid substructure found"

    # Look for the glycosidic linkage with at least one C-O-C connection
    # Custom pattern to check for likely glycosidic linkage
    connection_smarts = "O-[C@@H]"
    connection_pattern = Chem.MolFromSmarts(connection_smarts)

    if not mol.HasSubstructMatch(connection_pattern):
        return False, "No suitable glycosidic linkage found"
    
    return True, "Detected beta-D-glucuronic acid with a glycosidic linkage"