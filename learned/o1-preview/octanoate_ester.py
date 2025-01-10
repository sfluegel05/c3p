"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is an ester where the acid component is octanoic acid (caprylic acid),
    which has an 8-carbon unbranched, acyclic chain including the carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for octanoate ester
    # This pattern matches an ester group with an 8-carbon unbranched acyl chain
    octanoate_ester_smarts = '[#6;X3](=O)O[#6][#6][#6][#6][#6][#6][#6]'
    octanoate_ester_pattern = Chem.MolFromSmarts(octanoate_ester_smarts)

    if octanoate_ester_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Find matches in the molecule
    matches = mol.GetSubstructMatches(octanoate_ester_pattern)
    if matches:
        return True, "Contains octanoate ester group with 8-carbon unbranched acyl chain"
    else:
        return False, "No octanoate ester groups with unbranched 8-carbon acyl chain found"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'octanoate ester',
        'definition': 'Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).',
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}