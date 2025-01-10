"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.

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

    # SMARTS pattern for primary alcohol attached directly to an aromatic ring
    primary_alcohol_aromatic_pattern = Chem.MolFromSmarts("[#6]([CH2])[OH]")  # carbon with CH2 group attached to -OH

    # Find matches for the primary alcohol pattern in the molecule
    matches = mol.GetSubstructMatches(primary_alcohol_aromatic_pattern)

    for match in matches:
        # Check if the carbon is directly attached to an aromatic ring
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6:  # Check if aromatic neighbor is a carbon
                return True, "Contains a primary alcohol group directly attached to an aromatic ring"

    return False, "No primary alcohol group directly attached to an aromatic ring found"