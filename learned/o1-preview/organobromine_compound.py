"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carbon-bromine bond
    c_br_pattern = Chem.MolFromSmarts("[#6]-[#35]")
    if mol.HasSubstructMatch(c_br_pattern):
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "No carbon-bromine bond found"

__metadata__ = {   
    'chemical_class': {   
        'name': 'organobromine compound',
        'definition': 'A compound containing at least one carbon-bromine bond.',
    },
    'message': None,
    'success': True,
    'error': '',
}