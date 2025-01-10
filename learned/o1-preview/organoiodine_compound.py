"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: organoiodine compound
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound is defined as a compound containing at least one carbon-iodine (C-I) bond.

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

    # Define SMARTS pattern for carbon-iodine bond
    ci_pattern = Chem.MolFromSmarts("[#6]-[#53]")  # Carbon connected to iodine

    # Search for the pattern in the molecule
    if mol.HasSubstructMatch(ci_pattern):
        return True, "Contains at least one carbon-iodine bond"
    else:
        return False, "No carbon-iodine bonds found"