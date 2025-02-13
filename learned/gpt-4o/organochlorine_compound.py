"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains at least one carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carbon-chlorine bond
    c_cl_bond_pattern = Chem.MolFromSmarts("[#6]-[#17]")

    # Use pattern matching to find if there is a C-Cl bond
    if mol.HasSubstructMatch(c_cl_bond_pattern):
        return True, "Contains at least one carbon-chlorine bond"

    return False, "No carbon-chlorine bond found"