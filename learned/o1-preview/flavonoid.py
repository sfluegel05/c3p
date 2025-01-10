"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is any compound based on 1-benzopyran with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavonoid core SMARTS pattern
    # This pattern represents a 1-benzopyran ring fused to a benzene ring (A and C rings)
    # with an aryl group at position 2 (B ring)
    flavonoid_smarts = """
    [$([cH]1ccc2c(c1)occ2[$(c3ccccc3),$(c3ccc(cc3)O),$(c3ccc(cc3)OC),$(c3ccc(cc3)C)])]
    """

    # Remove whitespace and newlines from SMARTS
    flavonoid_smarts = ''.join(flavonoid_smarts.split())

    # Convert SMARTS to RDKit mol object
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smarts)
    if flavonoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Perform substructure search
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains flavonoid core structure with aryl substituent at position 2"
    else:
        return False, "Does not contain flavonoid core structure with aryl substituent at position 2"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47916',
                              'name': 'flavonoid',
                              'definition': 'Any member of the \'superclass\' '
                                            'flavonoids whose skeleton is based '
                                            'on 1-benzopyran with an aryl '
                                            'substituent at position 2. The '
                                            'term was originally restricted to '
                                            'natural products, but is now also '
                                            'used to describe semi-synthetic '
                                            'and fully synthetic compounds.',
                              'parents': []},
        'config': {},
        'message': None,
        'attempt': 1,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}