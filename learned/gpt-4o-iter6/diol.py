"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol should contain exactly two hydroxy groups positioned to be considered a diol, i.e., on adjacent or nearby carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a diol, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to detect hydroxy groups (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count the number of hydroxy groups
    num_hydroxy_groups = len(hydroxy_matches)

    # A diol should have exactly two hydroxy groups
    if num_hydroxy_groups != 2:
        return False, f"Contains {num_hydroxy_groups} hydroxy groups, need exactly 2"

    # Check connectivity of hydroxy groups to see if they are on adjacent or nearby carbons
    # Create a list to store the indices of the carbon atoms attached to hydroxy groups
    hydroxy_carbon_indices = []

    for match in hydroxy_matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                hydroxy_carbon_indices.append(neighbor.GetIdx())
    
    # Check connectivity of the hydroxy carbons
    # There should be a path or short distance within the molecule between the two hydroxy-bearing carbons
    if len(hydroxy_carbon_indices) == 2:
        path = Chem.rdmolops.GetShortestPath(mol, hydroxy_carbon_indices[0], hydroxy_carbon_indices[1])
        if len(path) <= 3:  # Allowing for diols like 1,2 or 1,3-diol
            return True, "Contains exactly two hydroxy groups forming a suitable diol"
    
    return False, "Hydroxy groups are not adequately positioned to form a diol"