"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is a compound where a carbonyl group is bonded to one hydrogen atom and to one R group (RC(=O)H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for aldehyde group: [CX3H1](=O)[#6]
    # General pattern where [CX3H1] is a carbon atom double bonded to oxygen and bonded to one hydrogen.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")

    # Search for the aldehyde pattern in the molecule
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Aldehyde functional group found"
    else:
        return False, "Aldehyde functional group not found"

# Test cases can be evaluated using this function to classify them as aldehydes.