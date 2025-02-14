"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by a sulfonate group (-S(=O)(=O)[O-])
    attached to an alkane chain (carbon chain).

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

    # Check for sulfonate group (-S(=O)(=O)[O-])
    sulfonate_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[O-]")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "No sulfonate group found"

    # Check for carbon directly attached to sulfonate sulfur
    carbon_sulfur_pattern = Chem.MolFromSmarts("[CX4][S](=[O])(=[O])[O-]")
    if not mol.HasSubstructMatch(carbon_sulfur_pattern):
        return False, "No carbon directly attached to sulfonate sulfur"

    return True, "Molecule contains a sulfonate group with at least one carbon attached to the sulfur atom"