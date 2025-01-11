"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is a compound containing at least one carbon-fluorine bond.

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

    # Define SMARTS pattern for carbon-fluorine bond
    cf_bond = Chem.MolFromSmarts("[#6]-[#9]")  # Carbon atom bonded to Fluorine atom
    if mol.HasSubstructMatch(cf_bond):
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "No carbon-fluorine bond found"

__metadata__ = {
    'chemical_class': {
        'id': '',
        'name': 'organofluorine compound',
        'definition': 'An organofluorine compound is a compound containing at least one carbon-fluorine bond.',
        'parents': []
    }
}