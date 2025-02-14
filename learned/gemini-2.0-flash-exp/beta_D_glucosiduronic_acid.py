"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is a glucuronic acid moiety with a glycosidic bond at the C1 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucuronic acid substructure with explicit stereochemistry and glycosidic bond at C1
    #  The * wildcard represents the O atom of the glycosidic bond.

    glucuronic_acid_smarts = "[C@H]1([C@@H]([C@H]([C@@H]([C@@H](O1)O[*])O)O)C(=O)O)"
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)


    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
         return False, "No beta-D-glucuronic acid substructure with glycosidic bond found"


    return True, "Contains beta-D-glucuronic acid with a glycosidic bond at the C1 position."