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

    # Define SMARTS pattern to find C-Cl bonds
    cl_pattern = Chem.MolFromSmarts("[#6]-[#17]")  # Carbon bonded to Chlorine

    # Search for the pattern in the molecule
    if mol.HasSubstructMatch(cl_pattern):
        return True, "Contains at least one carbon-chlorine bond"
    
    return False, "No carbon-chlorine bond found"