"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine compound
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
    
    # Look for carbon-fluorine bonds (C-F)
    cf_pattern = Chem.MolFromSmarts("[CX4,CX3][F]")
    if mol.HasSubstructMatch(cf_pattern):
        return True, "Contains at least one carbon-fluorine bond"

    return False, "No carbon-fluorine bonds found"


__metadata__ = {
    'chemical_class': {
        'id': 'Organofluorine compound',
        'name': 'organofluorine compound',
        'definition': 'An organofluorine compound is a compound containing at least one carbon-fluorine bond.'
    }
}