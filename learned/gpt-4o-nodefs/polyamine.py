"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is characterized by the presence of multiple amine groups connected across the structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extended amine pattern to include primary, secondary, tertiary, and cyclic amines
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;+1,0;!$([N]C=O)]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)

    # Minimum requirement is at least two amine-like groups, but also checking pattern distribution
    if len(amine_matches) < 2:
        return False, f"Found {len(amine_matches)} amine group(s), need at least 2"

    # Conditions for polyamine could include the dispersal across the molecule
    # Look for potential chain distribution across at least one long path
    # Check through adjacency to avoid recognizing non-chain configurations incorrectly
    # A heuristic sequence is necessary—identify more than 2 amines spaced across multiple bonds
    
    # Establish exceptions for missed connections (cyclic, special names, multi-center)
    amine_atom_idxs = {idx for match in amine_matches for idx in match}
    shortest_path_length = float('inf')
    
    # Identify atom path lengths to check dispersal span
    for i in amine_atom_idxs:
        for j in amine_atom_idxs:
            if i < j:  # Avoid redundant checks
                path_length = len(Chem.rdmolops.GetShortestPath(mol, i, j))
                if path_length > 4:  # Ensure a reasonable chain length
                    shortest_path_length = min(shortest_path_length, path_length)

    if shortest_path_length == float('inf') or shortest_path_length < 5:
        return False, "Amine groups not sufficiently spaced"
    
    return True, f"Contains multiple spaced amine groups with sufficient dispersion"

# The function is_polyamine can now be used to classify polyamines based on the SMILES input more accurately.