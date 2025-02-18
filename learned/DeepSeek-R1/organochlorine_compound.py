"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compounds (compounds with at least one carbon-chlorine bond)
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on the presence of at least one carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has a C-Cl bond, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a chlorine atom bonded to a carbon atom
    pattern = Chem.MolFromSmarts("[Cl]-[#6]")
    
    # Check for presence of the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains at least one carbon-chlorine bond"
    else:
        return False, "No carbon-chlorine bonds found"