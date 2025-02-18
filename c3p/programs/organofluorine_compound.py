"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:37143 organofluorine compound
An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined as a compound containing at least one carbon-fluorine bond.

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

    # Check for at least one carbon-fluorine bond using SMARTS pattern
    c_f_bond_pattern = Chem.MolFromSmarts("[C]-[F]")
    if mol.HasSubstructMatch(c_f_bond_pattern):
        return True, "Contains at least one carbon-fluorine bond"
    
    return False, "No carbon-fluorine bonds found"