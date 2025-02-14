"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:24842 indole alkaloid
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is an alkaloid containing an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general indole SMARTS pattern
    # This pattern matches indole rings, including substituted and fused ones
    indole_pattern = Chem.MolFromSmarts('n1c2cccc2cc1')  # Indole ring

    # Check for indole skeleton
    if mol.HasSubstructMatch(indole_pattern):
        # Check for presence of nitrogen atoms (since alkaloids contain nitrogen)
        atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        if 7 in atom_nums:
            return True, "Contains indole skeleton with nitrogen atoms (alkaloid)"
        else:
            return False, "Contains indole skeleton but no nitrogen atoms (not an alkaloid)"
    else:
        return False, "No indole skeleton found"