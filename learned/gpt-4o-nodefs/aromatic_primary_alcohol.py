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

    # Improved SMARTS pattern for primary alcohol attached directly to an aromatic ring
    primary_alcohol_aromatic_pattern = Chem.MolFromSmarts("[#6;H2][OH]")  # carbon with 2 hydrogens [H2] attached to -OH

    # Check if any aromatic carbon (aromatic atoms) is directly bonded to this primary alcohol group
    aromatic_carbons = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic() and atom.GetAtomicNum() == 6]
    primary_alcohol_matches = mol.GetSubstructMatches(primary_alcohol_aromatic_pattern)

    for carbon_atom in aromatic_carbons:
        carbon_neighbors = [neighbor.GetIdx() for neighbor in carbon_atom.GetNeighbors()]
        for match in primary_alcohol_matches:
            if match[0] in carbon_neighbors:
                return True, "Contains a primary alcohol group directly attached to an aromatic ring"

    return False, "No primary alcohol group directly attached to an aromatic ring found"