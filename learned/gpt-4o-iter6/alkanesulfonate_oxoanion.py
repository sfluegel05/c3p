"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by a sulfonate group directly attached
    to an alkane carbon atom.

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

    # A SMARTS pattern for an alkane carbon attached to a sulfonate group
    alkanesulfonate_pattern = Chem.MolFromSmarts("[CX4,CH2]-S(=O)(=O)[O-]")

    # Check if the sulfonate is directly attached to a carbon atom of type CH2 (alkane-like)
    if not mol.HasSubstructMatch(alkanesulfonate_pattern):
        return False, "Does not contain the characteristic alkanesulfonate sulfonate group directly attached to an alkane carbon"

    # Additional checks can be added here to exclude complex or aromatic attachments

    return True, "Contains a characteristic alkanesulfonate sulfonate group attached to an alkane carbon"

# Example test cases
print(is_alkanesulfonate_oxoanion("OC[C@H](O)CS([O-])(=O)=O"))  # Example input