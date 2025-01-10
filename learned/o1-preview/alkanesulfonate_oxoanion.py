"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""

from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a sulfonate group (-SO3-) attached to a carbon atom,
    where the carbon can be connected to any group (hydrogen, carbon chain, or other groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alkanesulfonate oxoanion
    # Pattern: carbon atom connected to sulfur atom that is doubly bonded to two oxygens and singly bonded to an oxygen with a negative charge
    pattern = Chem.MolFromSmarts('[C]-S(=O)(=O)[O-]')
    if mol.HasSubstructMatch(pattern):
        return True, "Contains alkanesulfonate oxoanion group"
    else:
        return False, "Does not contain alkanesulfonate oxoanion group"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'alkanesulfonate oxoanion',
        'definition': 'An alkanesulfonate in which the carbon at position 1 is attached to R, which can represent hydrogens, a carbon chain, or other groups.',
        'parents': []
    }
}