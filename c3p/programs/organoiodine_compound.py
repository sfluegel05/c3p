"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:? organoiodine compound
An organoiodine compound is defined as a compound containing at least one carbon-iodine bond.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound must contain at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organoiodine compound, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carbon-iodine bond.
    # [#6] is any carbon atom and [#53] is iodine. The '-' indicates a single bond.
    pattern = Chem.MolFromSmarts("[#6]-[#53]")
    
    # Check if the molecule has at least one carbon-iodine bond.
    if mol.HasSubstructMatch(pattern):
        return True, "Compound contains at least one carbon-iodine bond."
    else:
        return False, "Compound does not contain any carbon-iodine bonds."