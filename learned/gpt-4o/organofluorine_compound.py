"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound contains at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a C-F bond pattern that matches any carbon-fluorine connection
    cf_pattern = Chem.MolFromSmarts("[#6]-[#9]")
    
    # Check for C-F bonds
    if mol.HasSubstructMatch(cf_pattern):
        return True, "Contains at least one C-F bond"

    return False, "Does not contain any C-F bonds"