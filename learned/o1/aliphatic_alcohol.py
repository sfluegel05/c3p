"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: Aliphatic Alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound (no aromatic rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic atoms (molecule should not contain any)
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if aromatic_atoms:
        return False, "Molecule contains aromatic rings"

    # Define SMARTS pattern for hydroxyl group attached to carbon
    alcohol_pattern = Chem.MolFromSmarts('[CX4,CX3][OX2H]')

    # Find matches of the alcohol pattern
    matches = mol.GetSubstructMatches(alcohol_pattern)
    if not matches:
        return False, "No aliphatic hydroxyl groups found"

    # Ensure the hydroxyl groups are not attached to aromatic carbons
    for match in matches:
        carbon_idx = match[0]
        oxygen_idx = match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        if carbon_atom.GetIsAromatic():
            return False, "Hydroxyl group attached to aromatic carbon"

    return True, "Contains aliphatic alcohol group(s)"