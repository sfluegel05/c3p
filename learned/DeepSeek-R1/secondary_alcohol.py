"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohols (CHEBI:35678)
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has a hydroxyl group (-OH) attached to a carbon atom
    that is bonded to exactly two other carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all hydroxyl oxygen atoms (OX2H matches -OH groups)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    for match in hydroxyl_matches:
        oxygen_idx = match[0]
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        
        # Get the connected carbon (should be only one neighbor for -OH)
        neighbors = oxygen.GetNeighbors()
        if len(neighbors) != 1:
            continue  # Not a typical hydroxyl group
        
        carbon = neighbors[0]
        if carbon.GetAtomicNum() != 6:
            continue  # Hydroxyl not attached to carbon
        
        # Count carbon neighbors (should have exactly two other carbons)
        carbon_neighbors = [atom for atom in carbon.GetNeighbors() if atom.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 2:
            return True, "Hydroxyl attached to secondary carbon"

    return False, "No secondary alcohol group detected"