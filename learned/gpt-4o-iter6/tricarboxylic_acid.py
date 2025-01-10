"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined as having exactly three distinct and independent carboxylic acid groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tricarboxylic acid, False otherwise 
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid group using SMARTS
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Extract distinct carboxyl carbons
    distinct_carboxyl_carbons = set(match[0] for match in carboxylic_acid_matches)
    
    # Further inspect the groups to ensure they're individual:
    # We must ensure the co-existence of C(=O)O as an atom trajectory from the carboxyl carbon.
    valid_carboxyl_groups = []
    for carb_idx in distinct_carboxyl_carbons:
        atom = mol.GetAtomWithIdx(carb_idx)
        oxygen_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'O']
        
        # Ensure two oxygen atoms are directly bonded: one for =O, one for -OH
        if len(oxygen_neighbors) != 2:
            continue
        
        # Ensure no additional, non-hydrogen neighbors on carboxy(OH) beyond the initial carbon, oxygen, hydrogen
        for oxy_atom in oxygen_neighbors:
            more_neighbors_check = [
                neighbor for neighbor in oxy_atom.GetNeighbors() 
                if neighbor.GetSymbol() not in ['C', 'H', 'O']
            ]
            if len(more_neighbors_check) != 0:
                break
        else:
            valid_carboxyl_groups.append(carb_idx)

    # Count distinct and independent carboxylic acid groups
    if len(valid_carboxyl_groups) == 3:
        return True, "Contains exactly three distinct and independent carboxylic acid groups"
    else:
        return False, f"Contains {len(valid_carboxyl_groups)} distinct and independent carboxylic acid groups, expected exactly 3"