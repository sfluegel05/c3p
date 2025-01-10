"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by the presence of an -NH2 group where
    one hydrogen atom of ammonia (NH3) is replaced by a hydrocarbyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a primary amine group, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for primary amine group: nitrogen with two hydrogens and bonded to at least one carbon
    # SMARTS: [NH2][C], nitrogen atom bonded to exactly two hydrogens and bonded to a carbon
    primary_amine_pattern = Chem.MolFromSmarts("[NH2][CH2,CH,CH0,CX3,CX4]")
    if primary_amine_pattern is not None and mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Primary amine group detected"
    else:
        return False, "No primary amine group found"