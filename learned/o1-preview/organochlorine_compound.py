"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: An organochlorine compound is a compound containing at least one carbon-chlorine bond.
"""

from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound is a compound containing at least one carbon-chlorine bond.

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

    # Define the SMARTS pattern for a carbon-chlorine bond
    c_cl_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
    if mol.HasSubstructMatch(c_cl_pattern):
        return True, "Contains at least one carbon-chlorine bond"
    else:
        return False, "Does not contain any carbon-chlorine bond"