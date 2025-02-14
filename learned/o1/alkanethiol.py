"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thiol pattern: sulfur with one hydrogen attached to any carbon
    thiol_pattern = Chem.MolFromSmarts("[#16H1]-[#6]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)

    if not thiol_matches:
        return False, "Does not contain an -SH group attached to a carbon"

    # Check for disulfides or other sulfur-sulfur bonds
    disulfide_pattern = Chem.MolFromSmarts("[*]-[#16]-[#16]-[*]")
    if mol.HasSubstructMatch(disulfide_pattern):
        return False, "Contains disulfide bonds, not an alkanethiol"

    return True, "Contains a sulfanyl group (-SH) attached to an alkyl group"