"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carbon-iodine bond pattern: iodine directly bonded to carbon.
    carbon_iodine_pattern = Chem.MolFromSmarts("[I]-C")

    # Check if the molecule has a carbon-iodine bond
    if mol.HasSubstructMatch(carbon_iodine_pattern):
        return True, "Contains at least one carbon-iodine bond"
    else:
        return False, "Does not contain a carbon-iodine bond"