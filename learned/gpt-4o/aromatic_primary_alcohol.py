"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxyl (-OH) group attached to a carbon that is bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS query for aromatic primary alcohol
    # Aromatic ring carbon connected to -OH (alcohol) group
    aromatic_primary_alcohol_smarts = Chem.MolFromSmarts('[c][C](O)')

    # Check if the molecule matches the aromatic primary alcohol SMARTS pattern
    if mol.HasSubstructMatch(aromatic_primary_alcohol_smarts):
        return True, "Matches aromatic primary alcohol pattern"

    return False, "Does not match aromatic primary alcohol pattern"