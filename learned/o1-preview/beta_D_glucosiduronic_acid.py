"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is a molecule resulting from the formal condensation of any substance with beta-D-glucuronic acid to form a glycosidic bond.

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

    # Define the beta-D-glucuronic acid moiety SMARTS pattern
    # The anomeric carbon (C1) connected via O-glycosidic bond to any group ([!H])
    # Stereochemistry is specified using @@ and @ symbols
    # The carboxylic acid group at C6 can be protonated or deprotonated ([O-,OH])
    glucuronic_acid_smarts = '[C@@H]1([O][!H])[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C(=O)[O-,OH]'
    
    # Create the glucuronic acid pattern molecule
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)
    if glucuronic_acid_pattern is None:
        return False, "Failed to create beta-D-glucuronic acid pattern"

    # Check for the beta-D-glucuronic acid substructure
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    return True, "Contains beta-D-glucuronic acid moiety linked via glycosidic bond"