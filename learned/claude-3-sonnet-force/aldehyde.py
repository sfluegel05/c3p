"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:36975 aldehyde
An aldehyde is a compound RC(=O)H, in which a carbonyl group is bonded to one hydrogen atom and to one R group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns to identify aldehyde groups
    aldehyde_patterns = [
        "[H]C=O",  # Aldehyde group with explicit hydrogen
        "[H]/C=O", # Aldehyde group with explicit hydrogen (forward slash)
        "[H]\\C=O", # Aldehyde group with explicit hydrogen (backward slash)
        "C(=O)[H]", # Aldehyde group with explicit hydrogen
    ]

    # Check if the molecule contains any of the aldehyde patterns
    for pattern in aldehyde_patterns:
        aldehyde_query = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(aldehyde_query)

        # Check if any matches are valid aldehyde groups
        for match in matches:
            if is_valid_aldehyde_group(mol, match):
                return True, "Molecule contains a valid aldehyde group"

    return False, "No valid aldehyde group found in the molecule"

def is_valid_aldehyde_group(mol, match):
    """
    Checks if a given match in the molecule corresponds to a valid aldehyde group.

    Args:
        mol (Mol): RDKit molecule object
        match (tuple): Tuple of atom indices representing the match

    Returns:
        bool: True if the match is a valid aldehyde group, False otherwise
    """
    carbonyl_idx = match[0]  # Assume the first atom in the match is the carbonyl carbon

    # Check if the carbonyl carbon has exactly one non-hydrogen neighbor
    non_hydrogen_neighbors = [
        nbr_idx for nbr_idx in mol.GetAtomWithIdx(carbonyl_idx).GetNeighbors()
        if mol.GetAtomWithIdx(nbr_idx).GetSymbol() != "H"
    ]

    if len(non_hydrogen_neighbors) != 1:
        return False

    # Check if the non-hydrogen neighbor is not part of a ring system
    neighbor_idx = non_hydrogen_neighbors[0]
    if mol.GetAtomWithIdx(neighbor_idx).IsInRing():
        return False

    # Additional checks for ring systems, functional groups, etc. can be added here

    return True