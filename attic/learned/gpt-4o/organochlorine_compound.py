"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound must contain at least one carbon-chlorine bond.

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

    # Look for carbon-chlorine bond pattern
    chlorocarbon_pattern = Chem.MolFromSmarts("[#6]-[#17]")  # C-Cl pattern
    if mol.HasSubstructMatch(chlorocarbon_pattern):
        return True, "Contains at least one carbon-chlorine bond"
        
    return False, "No carbon-chlorine bond found"