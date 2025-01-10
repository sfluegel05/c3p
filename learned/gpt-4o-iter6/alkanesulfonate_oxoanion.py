"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by a sulfonate group directly attached
    to an alkyl carbon atom, where R can be hydrogen, a carbon chain, or other groups.

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

    # A SMARTS pattern for carbon (R) attached to a sulfonate group (C-SO3)
    alkanesulfonate_pattern = Chem.MolFromSmarts("[CX3,CX4][SX3](=O)(=O)[O-]")

    # Match the pattern in the molecule
    if not mol.HasSubstructMatch(alkanesulfonate_pattern):
        return False, "Does not contain the characteristic alkanesulfonate sulfonate group attached to an alkyl carbon"

    # Further checks can be added here if required to differentiate more complex aromatic systems
    # E.g., exclude if part of an aromatic ring
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S' and atom.GetNeighbors():
            is_aromatic = any(neighbor.GetIsAromatic() for neighbor in atom.GetNeighbors())
            if is_aromatic:
                return False, "Part of a complex aromatic system, not a typical alkanesulfonate"

    return True, "Contains a characteristic alkanesulfonate sulfonate group attached to an alkyl carbon"

# Example test cases
print(is_alkanesulfonate_oxoanion("OC[C@H](O)CS([O-])(=O)=O"))  # Example input