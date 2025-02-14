"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:16445 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as an oxime derived from an aldehyde, having the functional group R-CH=N-OH,
    where the carbon double-bonded to nitrogen is attached to exactly one hydrogen, distinguishing it from
    ketoximes where the carbon is attached to two carbon groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldoxime functional group (R-CH=N-OH)
    aldoxime_pattern = Chem.MolFromSmarts('[CX3H1]=N[OH]')
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "No aldoxime functional group found"

    return True, "Contains aldoxime functional group (-CH=N-OH)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16445',
        'name': 'aldoxime',
        'definition': 'Oximes of aldehydes RCH=NOH.',
        'parents': []
    },
}