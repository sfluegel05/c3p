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
    An organoiodine compound contains at least one carbon-iodine (C-I) bond.

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
    ci_bond_pattern = Chem.MolFromSmarts("[#6]-[#53]")
    if ci_bond_pattern is None:
        return None, "Error creating SMARTS pattern for C-I bond"

    # Search for C-I bonds
    matches = mol.GetSubstructMatches(ci_bond_pattern)
    if matches:
        return True, f"Found {len(matches)} carbon-iodine bond(s)"
    else:
        return False, "No carbon-iodine bonds found"