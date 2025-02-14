"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is a glycerol backbone with three fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester pattern: carbonyl carbon(=O) - ester oxygen - glycerol carbon
    ester_pattern = Chem.MolFromSmarts('[#6](=O)O[#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Collect the carbons attached to ester oxygens (glycerol carbons)
    esterified_carbons = set()
    for match in ester_matches:
        ester_o_idx = match[1]      # Index of the ester oxygen atom
        glycerol_c_idx = match[2]   # Index of the glycerol carbon attached to ester oxygen
        esterified_carbons.add(glycerol_c_idx)

    # Check if there are exactly 3 esterified carbons
    if len(esterified_carbons) != 3:
        return False, f"Found {len(esterified_carbons)} ester groups, need exactly 3"

    # Convert the set to a list for indexing
    glycerol_carbons = list(esterified_carbons)

    # Check if the esterified carbons form a connected chain (C-C-C)
    # Compute the shortest paths between the carbons
    path_lengths = []
    for i in range(3):
        for j in range(i+1, 3):
            path = Chem.rdmolops.GetShortestPath(mol, glycerol_carbons[i], glycerol_carbons[j])
            path_lengths.append(len(path) - 1)  # Subtract 1 to get number of bonds

    # For a linear chain C-C-C, the path lengths should be [1, 1, 2]
    if sorted(path_lengths) != [1, 1, 2]:
        return False, "Esterified carbons do not form a glycerol backbone (C-C-C chain)"

    # If all checks pass, the molecule is a triglyceride
    return True, "Contains esterified glycerol backbone with three fatty acid chains"